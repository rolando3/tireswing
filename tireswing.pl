#! /usr/bin/perl

use strict;
use Getopt::Long;
use Digest::MD5 qw(md5 md5_hex);
use POSIX qw(floor ceil);

#these are the three parameters:
# 1. Who are my focal individuals
# 2. What't the threshold; i.e. how many links are required?
# 3. How many times do we want to iterate?

#these are parameters
my ( @foci, $foci ); #focal nodes. Should be names as expected in input. Must match exactly
my ( @exclude, $exclude ) = ( undef, {} ); #individuals to exclude from the map. If individual
                                           #is in @foci then this just prevents printing. if
                                           #not in foci then this individual is ignored when
                                           #walking the graph
my $distance=1; #Maximum iterations
my $keep=1; #keep the focal individuals
my $keepnew=1; #on the last iteration, find links between newly discovered individuals
my $verbose; #if true sends warnings of progress to stderr
my $chr; #limit to a single chromosome
my ( $range, @range ); #limit to a range on a given chromosome (e.g., 29-35)
my $outputmode = "s"; #Type of output. s = one row per segment, p = one row per person/person match
my $max; #maximum number of vertices to output - not exactly what we do
my $hide; #if true then hashes names to protect the innocent
my $thresholdoverride;
my $exclusive;

my $help;

sub usage
{
    die "perl afscan.pl -f focal_individual [-x exclude_individual] [-d distance] [-o s|p] [-keep|nokeep] [-new|nonew] [-max maxmatches] [-chr chromosome [-range range]] [-hide] filelist\n";
}

usage() if ( @ARGV < 1 or
    ! GetOptions(
        'focal|f=s' => \@foci,
        'exclude|x:s' => \@exclude,
        'distance|d:i' => \$distance,
        'output|o:s' => \$outputmode,
        'keep|k!' => \$keep,
        'new|n!' => \$keepnew,
        'max|m:i' => \$max,
        'verbose|v+' => \$verbose,
        'chr|c:s' => \$chr,
        'range:s' => \$range,
        'hide|h' => \$hide,
        'exclusive|e!' => \$exclusive,
        'thresholdoverride|to:i' => \$thresholdoverride,
        'help|?' => \$help
    ) or $help );

unless ( ( $outputmode == "s" ) or ( $outputmode == "p" ) ) {
    die "Invalid output mode. Valid modes are 's' (segments) and 'p' (people).\n"
}

if ( $range ) {
    die "Range must be accompanied by chromosome specification.\n" unless $chr;
    foreach ( $range ) { @range = split "-"; }
    die "Range ($range) must be specified as start-end.\n" unless ( @range == 2 );
    warn " $range[0] - $range[1]\n";
}

map { $exclude->{$_} = 1; } @exclude;
map { $foci->{$_} = 1; } @foci;

my $rawdata = {};
my $newmatches = {};
my $estdmatches = {};
my $outputdata = {};

sub read_data
{
    my $uid;
    while (<>)
    {
        chomp;
        $uid ++;
        my @l = split ",";
        map { s/\"//g; } @l;

        my $p1 = mask_name($l[0]);
        my $p2 = mask_name($l[1]);

        my $link = {
            p1 => $p1,
            p2 => $p2,
            chr => $l[10],
            start => $l[11],
            end => $l[12],
            cM => $l[14],
            id => $uid
        };

        if ( $chr )
        {
            next if ( $link->{chr} ne $chr );
            if ( $range )
            {
                next if ( ( $link->{start} > $range[1] ) or ( $link->{end} < $range[0] ) );
            }
        }

        #we are making a hash of arrays - each a list of matches for the name in p1
        $rawdata->{$p1} = {} unless ( $rawdata->{$p1} );
        $rawdata->{$p1}->{$p2} = {} unless ( $rawdata->{$p1}->{$p2} );

        $rawdata->{$p2} = {} unless $rawdata->{$p2};
        $rawdata->{$p2}->{$p1} = {} unless ( $rawdata->{$p2}->{$p1} );
        
        my $loc = "$link->{chr}.$link->{start}";
        
        #add a match to p1's array of matches
        $rawdata->{$p1}->{$p2}->{$loc} = $link;
        $rawdata->{$p2}->{$p1}->{$loc} = $link;

    }

    warn "Read $uid matches...\n" if $verbose;
}

sub add_foci
{
    foreach ( @foci )
    {
        if ( $rawdata->{$_} )
        {
            $estdmatches->{$_} = 1;
        }
        else
        {
            warn "$_ not found.\n";
        }
    }

    my @f = keys(%$estdmatches);
    die "No valid foci given. Check spelling and underscores/spaces.\n" if ( @f eq 0 );

    warn "Added ".@f." focal match(es)...\n" if $verbose;

    map { add_to_output ( $_, find_links ( $_, $estdmatches ) ); } @f;

    $newmatches = find_candidate_links ( @f )
}

sub add_to_output
{
    my ( $nm, @t ) = @_;
    #if so then go through his matches.
    foreach ( @t )
    {
        die "what is going on here?\n" if ($nm eq $_);
        if ( $nm gt $_ ) {
            $outputdata->{$_} = {} unless $outputdata->{$_};
            $outputdata->{$_}->{$nm} = 1;
        } else {
            $outputdata->{$nm} = {} unless $outputdata->{$nm};
            $outputdata->{$nm}->{$_} = 1;
        }
    }

}

sub find_links
{
    #look through the target list and return a list of people who match this guy.
    my ( $source, $targetlist ) = @_;
    my @t;

    foreach ( keys %$targetlist )
    {
        next if ( $_ eq $source );
        push ( @t, $_ ) if ( $rawdata->{$source}->{$_} );
    }

#    print uc($source)." @t\n";
    return @t;
}

sub find_candidate_links
{
    my $result = {};
    while ( my $ind = shift )
    {
        foreach ( find_links ( $ind, $rawdata->{$ind} ) )
        {
            $result->{$_} = 1 unless ( $estdmatches->{$_} or $exclude->{$_} );
        }
    }

    return $result;
}

sub find_links_loop
{
    my $loopcount = shift;
    my $islastloop = ( $loopcount == $distance );
    my $newlyestablished = {};
    my @newmatchlist = keys %$newmatches;

    my $threshold = get_threshold($loopcount);
    warn "Scanning ".@newmatchlist." new matches (threshold: $threshold)...\n" if $verbose;

    #go through each new match.
    foreach ( @newmatchlist )
    {

        my $t = {};
        my $nm = $_;

        map { $t->{$_} = 1; } find_links ( $nm, $estdmatches );

        #does this new match meet the criteria for being saved?
        if ( keys(%$t) >= $threshold )
        {
#            warn "$nm has ".keys(%$t)." matches\n";
            add_to_output ( $nm, keys(%$t) );
            $newlyestablished->{$nm} = 1;
        }

    }

    my $newnewmatches = find_candidate_links(@newmatchlist);

    #move people along
    foreach ( keys ( %$newlyestablished ) )
    {
        $estdmatches -> {$_} = 1;
        delete $newmatches->{$_};
        add_to_output ( $_, find_links ( $_, $newlyestablished ) ) if ( $islastloop and $keepnew );
    }

    warn "Found ".keys(%$newlyestablished)." new matches (total: ".keys(%$estdmatches).")\n" if $verbose;

    foreach ( keys(%$newnewmatches) )
    {
        $newmatches->{$_} = 1;
    }

}

sub _name
{
    my $n = shift;
    return ( $hide and not $foci->{$n} ) ? substr(md5_hex($n),0,8) : $n;
}

sub print_results
{
    #print the header
    print join(",", ( $outputmode eq "s" ? qw/p1 p2 chr start end cM/ : qw/p1 p2 count cM/ ) )."\n";
    foreach ( sort ( keys ( %$outputdata ) ) )
    {
        next if ( ! $keep and $foci->{$_} );
        next if $exclude->{$_};
        my $m1 = $_;
        foreach ( sort ( keys ( %{$outputdata->{$m1}} ) ) )
        {
            next if ( (! $keep) and $foci->{$_} );
            next if $exclude->{$_};
            my $m2 = $_;
            my $segcount = 0;
            my $cMtotal = 0;

            foreach ( keys %{$rawdata->{$m1}->{$m2}} )
            {
                my $link = $rawdata->{$m1}->{$m2}->{$_};

                if ( $outputmode eq "s" )
                {
                    print join(",",(_name($m1),_name($m2),$link->{chr},$link->{start},$link->{end},$link->{cM}))."\n";
                }
                else
                {
                    $segcount ++;
                    $cMtotal += $link->{cM};
                }
            }

            print join(",",(_name($m1),_name($m2),$segcount,$cMtotal))."\n" if ( $outputmode eq "p" );
        }
    }
}

sub mask_name
{

    my $name = shift;

    foreach ( $name ) {
        s/ /_/g;
        s/-/_/g;
        s/[0-9]/_/g;
        s/'/_/g;
        s/\./_/g;
    }

    return $name;
}

sub get_threshold
{
    return $thresholdoverride if $thresholdoverride;

    my $loop = shift;
    my $matches = keys(%$estdmatches);

    if ( $loop <= 1 )
    { 
        return ( $matches <= 2 ) ? $matches : ceil($matches/2 + $exclusive) 
    } 
    else
    {
	return ( $matches <= 2 ) ? $matches : floor( log($matches) / ( $exclusive ? log(2) : 1 ) );
    }
}

sub main
{
    #match
    read_data();
    #gets our focal nodes and establishes connections. then
    add_foci();
    #work our way out
    for ( my $i = 1; $i <= $distance; $i++ )
    {
        find_links_loop ( $i );
        if ( $max and ( keys(%$estdmatches) > $max ) )
        {
            warn "Maximum node count ($max) exceeded (".keys(%$estdmatches).").\n";
            last;
        }
    }
    #print all results
    print_results();
}

main();
