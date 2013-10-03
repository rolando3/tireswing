#! /bin/bash

header="FileMatch","MatchName","MaternalGrandmotherBirthCountry","MaternalGrandfatherBirthCountry","PaternalGrandmotherBirthCountry","PaternalGrandfatherBirthCountry","MaternalGrandmotherDeclaredAshkenazi","MaternalGrandfatherDeclaredAshkenazi","PaternalGrandmotherDeclaredAshkenazi","PaternalGrandfatherDeclaredAshkenazi","Chromosome","SegmentStartInMegaBasePairs","SegmentEndInMegaBasePairs","SegmentLengthInMegaBasePairs","SegmentLengthInCentiMorgans"

echo $header
find $1 -name ancestry_finder*csv | xargs egrep -v "Anonymous|MatchName" | sed "s/^.*ancestry_finder_/\"/; s/_20[0-9][0-9][0-9][0-9][0-9][0-9].csv:/\",/"
