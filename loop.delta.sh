#!/bin/tcsh

# 1 D = 0.393430307 e*a_0
# 2.5417462310548435


set parent=`pwd`
set chgdiff=~/scripts/vtstscripts/chgdiff.pl

foreach r ( ???_DOS )

 echo $r
 cd $r

 cd z

 pwd

# $chgdiff CHGCAR chg/sur/CHGCAR
# mv CHGCAR_diff CHGCAR.noslab
# $chgdiff CHGCAR.noslab chg/mol/CHGCAR
# mv CHGCAR_diff dNCAR
# rm CHGCAR.noslab

 rm -f planechg.in
 ln -s ../../dncar.in planechg.in
 /global/home/users/itamblyn/nano1/butanediamine_junction/chgplane.x
 rm -f planechg.in

# echo "READ NOTE THIS FILE ABOUT DELETING PART OF planechg.dat"

 awk '\!/#/{i+=((($2 - 14.840528)/0.529177)*$3*2.5417462310548435); print $2, i}' trimmed.dat # i deleted the parts of the file that were not the molecule!!!!!! 

 cd ../..

end

