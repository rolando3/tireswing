#! /usr/bin/perl

use strict;
use Getopt::Long;
use Digest::MD5 qw(md5 md5_hex);
use POSIX qw(floor ceil);

#these are the three parameters:
# 1. Who are my focal individuals
# 2. What't the threshold; i.e. how many links are required?
# 3. How many times do we want to iterate?

my $rawdata = {};
my $estdmatches = {};
my $outputdata = {};
my $outputcounts = {};
my $affiles = {};
my $uid;

#these are parameters
my ( @foci, $foci ); #focal nodes. Should be names as expected in input. Must match exactly
my ( @exclude, $exclude ); #individuals to exclude from the map. If individual
                           #is in @foci then this just prevents printing. if
                           #not in foci then this individual is ignored when
                           #walking the graph
my $distance=1; #Maximum iterations
my $keepfocal=1; #keep the focal individuals
my $keepnew=1; #on the last iteration, find links between newly discovered individuals
my $header = 1;
my $verbose; #if true sends warnings of progress to stderr
my $chr; #limit to a single chromosome
my ( $range, @range ); #limit to a range on a given chromosome (e.g., 29-35)
my $outputmode = "s"; #Type of output. s = one row per segment, p = one row per person/person match
my $nodelimit; #maximum number of vertices to output - not exactly what we do
my $hide; #if true then hashes names to protect the innocent
my $thresholdoverride;
my $exclusive;
my $minsegmentlength;
my $multimatch;

my $help;

sub usage
{
    die "perl $0 -f focal_individual [-x exclude_individual] [-d distance] [-o s|p] [-multimatch] [-limit limit] [-verbose] [-chr chromosome [-range range]] [-hide] [-exclusive] [-thresholdoverride num] -- filelist\n";
}
#[-keep|nokeep] [-new|nonew]

sub main
{
    #get parameters from command line (or wherever?)
    get_params();

    my $starttime = time();
    my $newmatches = {};

    #match
    read_data(@ARGV);

    warn "Time to read data: ".(time - $starttime)."\n" if $verbose;

    #gets our focal nodes and establishes connections between them
    $newmatches = add_foci();
    #work our way out
    for ( my $i = 1; $i <= $distance; $i++ )
    {
        $newmatches = extend_distance ( $i, $newmatches );
        if ( $nodelimit and ( keys(%$estdmatches) > $nodelimit ) )
        {
            warn "Maximum node count ($nodelimit) exceeded (".keys(%$estdmatches).").\n";
            last;
        }
    }

    #print all results
    print_results();

    warn "Time elapsed: ".(time - $starttime)."\n" if $verbose;

}

sub get_params
{
    usage() if ( @ARGV < 1 or
        ! GetOptions(
            'focal|f=s' => \@foci,
            'exclude|x:s' => \@exclude,
            'distance|d:i' => \$distance,
            'output|o:s' => \$outputmode,
            'limit|l:i' => \$nodelimit,
            'verbose|v+' => \$verbose,
            'min-segment-length|msl:f' => \$minsegmentlength,
            'multimatch|mm' => \$multimatch,
            'chr|c:s' => \$chr,
            'range:s' => \$range,
            'hide|h' => \$hide,
            'exclusive|e!' => \$exclusive,
            'thresholdoverride|to:i' => \$thresholdoverride,
            'keepnew!' => \$keepnew,
            'keepfocal!' => \$keepfocal,
            'header!' => \$header,
            'help|?' => \$help
        ) or $help );

    unless ( ( $outputmode == "s" ) or ( $outputmode == "p" ) ) {
        die "Invalid output mode. Valid modes are 's' (segments) and 'p' (people).\n"
    }

    if ( $range ) {
        die "Range must be accompanied by chromosome specification.\n" unless $chr;
        foreach ( $range ) { @range = split "-"; }
        die "Range ($range) must be specified as start-end.\n" unless ( @range == 2 );
    }

    map { $exclude->{$_} = 1; } @exclude;
    map { $foci->{$_} = 1; } @foci;

    if ( $verbose )
    {
        warn "Parameters:\n";
        warn "  Foci: ".join(", ",@foci)."\n";
        warn "  Exclude: ".join(", ",@exclude)."\n" if @exclude;
        warn "  Distance: $distance\n";
        warn "  Output mode: $outputmode\n";
        warn "  Limit: ".($nodelimit || "N/A")."\n";
        warn "  Minimum segment length: ".($minsegmentlength||"N/A")."\n";
        warn "  Chromosome: ".($chr || "ALL")."\n";
        warn "  Chromosome Range: ".($range || "N/A")."\n";
        warn "  Limit output to nodes with multiple matches: $multimatch\n";
        warn "  Hash non-focal individuals' names: ".($hide||"N/A")."\n";
        warn "  Exclusive thresholds: ".($exclusive||"N/A")."\n";
        warn "  Threshold overrides: ".($thresholdoverride||"N/A")."\n";
    }

}

sub read_line
#reads a single line from a 23andMe ancestry finder output file and adds it to our dataset
{
    #two parameters: the owner of the file and the line itself
    my $u1 = _23andme_name_mask(shift);

    #break the line into an array
    chomp;
    my @l = split ",";
    map { s/\"//g; } @l;

    #make sure the name in the file has the same format as the owner name (underscores)
    my $u2 = _23andme_name_mask($l[0]);
    $uid ++;

    #build a link
    my $link = {
        u1 => $u1,
        u2 => $u2,
        chr => $l[9],
        start => $l[10],
        end => $l[11],
        cM => $l[12],
        id => $uid
    };

    #discard this link if check_link says it's not worthy
    next unless check_link ( $link );

    #we are making a hash of arrays - each a list of matches for the name in u1
    $rawdata->{$u1} = {} unless ( $rawdata->{$u1} );
    $rawdata->{$u1}->{$u2} = {} unless ( $rawdata->{$u1}->{$u2} );

    #do the same for u2. this is because both users may have files and may match each other
    #we are conservative in what we store. it might be possible to do this more judiciously
    $rawdata->{$u2} = {} unless $rawdata->{$u2};
    $rawdata->{$u2}->{$u1} = {} unless ( $rawdata->{$u2}->{$u1} );

    my $loc = "$link->{chr}.$link->{start}";

    #add a match to u1's array of matches
    $rawdata->{$u1}->{$u2}->{$loc} = $link;
    $rawdata->{$u2}->{$u1}->{$loc} = $link;
}

sub read_data
{
    warn "Opening @_ files..." if ( $verbose > 1 );
    #given a list of files, read the lines from each
    while ( $_ = shift )
    {
        #we extract the user name from each file name, which follow a particular format
        my $fn = $_;
        s/.*ancestry_finder_//;
        s/_2[0-9]{7}.csv//;
        my $name = $_;
        $affiles->{$name} = 1;

        warn "Opening $fn...\n" if ( $verbose > 1 );
        open(INPUT,$fn);

        #read each line in each file
        while (<INPUT>)
        {
            read_line $name, $_ unless m/Anonymous|MatchName/;
        }

        close INPUT;
    }
    warn "Read $uid matches...\n" if $verbose;
}

sub add_foci
{
    my @missing;
    foreach ( @foci )
    {
        if ( $rawdata->{$_} )
        {
            $estdmatches->{$_} = 1;
        }
        else
        {
            push @missing, $_;
        }
    }

    die "Invalid foci: ".join(",",@missing).". Check spelling and underscores/spaces.\n" if ( @missing );

    my @f = keys(%$estdmatches);

    warn "Added ".@f." focal match(es)...\n" if $verbose;

    map { add_to_output ( $_, find_links ( $_, $estdmatches ) ); } @f;

    return find_candidate_links ( @f )
}

sub add_to_output
{
    my ( $nm, @t ) = @_;
    #if so then go through his matches.
    foreach ( @t )
    {
        die "what is going on here?\n" if ($nm eq $_);
        $outputcounts->{$nm} ++;
        $outputcounts->{$_} ++;
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

    return @t;
}

sub find_candidate_links
{
    #given a list of names, find all the links that may match that list
    my $result = {};
    foreach ( @_ )
    {
        foreach ( find_links ( $_, $rawdata->{$_} ) )
        {
            $result->{$_} = 1 unless ( $estdmatches->{$_} or $exclude->{$_} );
        }
    }

    return $result;
}

sub check_link
{
    my $link = shift;
    #move along when user has specified a specific chromosome/range and this one doesn't match
    if ( $chr )
    {
        return 0 if ( $link->{chr} ne $chr );
        if ( $range )
        {
            return 0 if ( ( $link->{start} > $range[1] ) or ( $link->{end} < $range[0] ) );
        }
    }

    if ( $minsegmentlength )
    {
        return 0 if ( $link->{cM} < $minsegmentlength );
    }

    return 1;
}


sub extend_distance
{
    my ( $loopcount, $newmatches ) = @_;
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

        #for each new match we are going to find out which established matches they match
        map { $t->{$_} = 1; } find_links ( $nm, $estdmatches );

        #does this new match meet the criteria for being saved?
        if ( keys(%$t) >= $threshold )
        {
            #add it to the output list and mark it as "newly established"
            add_to_output ( $nm, keys(%$t) );
            $newlyestablished->{$nm} = 1;
        }

    }

    #move people along
    foreach ( keys ( %$newlyestablished ) )
    {
        #this is now an established match. will be included in the output unless we're not not keeping new ones
        $estdmatches -> {$_} = 1;
        #remove them from our list of new matches.
        delete $newmatches->{$_};
        #if this is the last loop and we're keeping new -- find their links among this group.
        add_to_output ( $_, find_links ( $_, $newlyestablished ) ) if ( $islastloop and $keepnew );
    }

    warn "Found ".keys(%$newlyestablished)." new matches (total: ".keys(%$estdmatches).")\n" if $verbose;

    #now it's time to put together our list of candidate links that match our current new match list
    map { $newmatches->{$_} = 1; } keys(%{find_candidate_links(@newmatchlist)}) unless $islastloop;

    return $newmatches;

}

sub print_results
{

    map { warn "Excluding $_ \n" if $outputcounts->{$_} == 1; } sort(keys(%$outputcounts)) if ( $verbose and $multimatch );

    #print the header
    print join(",", ( $outputmode eq "s" ? qw/u1 u2 chr start end cM u1focal u1data u2focal u2data/ : qw/u1 u2 count cM u1focal u1data u2focal u2data/ ) )."\n" if ( $header );
    foreach ( sort ( keys ( %$outputdata ) ) )
    {
        next if ( ! $keepfocal and $foci->{$_} );
        next if $exclude->{$_};
        next if ( $multimatch and $outputcounts->{$_} == 1 );

        my $u1 = $_;
        foreach ( sort ( keys ( %{$outputdata->{$u1}} ) ) )
        {
            next if ( ! $keepfocal and $foci->{$_} );
            next if $exclude->{$_};
            next if ( $multimatch and $outputcounts->{$_} == 1 );

            my $u2 = $_;
            my $segcount = 0;
            my $cMtotal = 0;

            my $on1 = get_output_name($u1);
            my $on2 = get_output_name($u2);

            foreach ( keys %{$rawdata->{$u1}->{$u2}} )
            {
                my $link = $rawdata->{$u1}->{$u2}->{$_};

                if ( $outputmode eq "s" )
                {
                    print join(",",($on1,$on2,$link->{chr},$link->{start},$link->{end},$link->{cM},is_focal($u1),has_data($u1),is_focal($u2),has_data($u2)))."\n";
                }
                else {
                    $segcount ++;
                    $cMtotal += $link->{cM};
                }
            }

            print join(",",$on1,$on2,$segcount,$cMtotal,is_focal($u1),has_data($u1),is_focal($u2),has_data($u2))."\n" if ( $outputmode eq "p" );
        }
    }
}

sub is_focal
{
    my $n = shift;
    return $foci->{$n} ? "1" : "0";
}

sub has_data
{
    my $n = shift;
    return $affiles->{$n} ? 1 : 0;
}

sub _23andme_name_mask
{
    foreach ( $_[0] ) {
        s/[ -'\.]/_/g;
        s/ñ/n/g;
        s/ó/o/g;
    }

    return $_[0];
    #        s/Ñ/N/g;
    #        s/á/a/g;
    #        s/ä/a/g;
    #        s/â/a/g;
    #        s/à/a/g;
    #        s/ã/a/g;
    #        s/å/a/g;
    #        s/é/e/g;
    #        s/ë/e/g;
    #        s/ê/e/g;
    #        s/è/e/g;
    #        s/í/i/g;
    #        s/ï/i/g;
    #        s/î/i/g;
    #        s/ì/i/g;
    #        s/ö/o/g;
    #        s/ô/o/g;
    #        s/ò/o/g;
    #        s/ø/o/g;
    #        s/[úüûù]/u/g;
    #        s/[ÁÄÂÀÃÅ]/A/g;
    #        s/[ÉËÊÈ]/E/g;
    #        s/[ÍÏÍÌ]/I/g;
    #        s/[ÓÖÔÒØ]/O/g;
    #        s/[ÚÜÛÙ]/U/g;
}

sub get_output_name
{
    my $o = shift;
    my $_ = $o;
    s/_/ /g;
    s/\W/ /g;              # convert others to spaces
    s/^\s+//;              # remove leading spaces
    s/\s+$//;              # remove trailing spaces

    my $rv = ( $hide and not $foci->{$o} ) ? substr(md5_hex($_),0,8) : $_;
    return cc($rv);
#    $_ = ( $foci->{$o} ) ? uc($_) : cc($_);
#    $_ = "* $_ *" if $affiles->{$o};
    return $_;
}

sub cc
{
    $_ = lc(shift);
    s/\s+(.)/ \u$1/g;       # remove other spaces and upcase next letter
    s/^(.)/\u$1/;          # downcase first letter
    return $_;
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
        my $t = ( $matches <= 2 ) ? $matches : floor( log($matches) / ( $exclusive ? log(2) : 1 ) );
        return ( $t <= 1 ) ? 2 : $t;
    }

    return 2;
}

main();
