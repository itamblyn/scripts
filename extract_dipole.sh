#!/bin/tcsh

# 1 D = 0.393430307 e*a_0
# 2.5417462310548435

rm -f planechg.in
ln -s ~/scripts/chg.in planechg.in
/global/home/users/itamblyn/bin/planechg.x
rm planechg.in planechg.out
mv planechg.dat chg.dat

rm -f planechg.in
ln -s ~/scripts/pot.in planechg.in
/global/home/users/itamblyn/bin/planechg.x
rm planechg.in planechg.out
mv planechg.dat pot.dat

module unload gnuplot
module load gnuplot/4.2.5-gcc-gd

set cmd="gnuplot.scr"
cat > $cmd << END
set style data lines
set xlabel "Position [A]"
set yrange [*:*]
set terminal png
set output 'density_and_potential.png'
plot "chg.dat" t "Charge density", "pot.dat" u 1:(\$2*100) t "Potential"
quit
END

gnuplot $cmd

\rm $cmd

if ( -d chg ) then

  if( ! -e trimmed.dat ) then

   set parent=`pwd`
   set chgdiff=~/scripts/vtstscripts/chgdiff.pl
   $chgdiff CHGCAR chg/sur/CHGCAR
   mv CHGCAR_diff CHGCAR.noslab
   $chgdiff CHGCAR.noslab chg/mol/CHGCAR
   mv CHGCAR_diff dNCAR
   rm CHGCAR.noslab
   rm -f planechg.in
   ln -s ~/scripts/dncar.in planechg.in
   /global/home/users/itamblyn/nano1/butanediamine_junction/chgplane.x
   rm -f planechg.in
   echo
   echo "Okay, now you have to trim the file planechg.dat, and call it trimmed.dat"
 endif

  if( -e trimmed.dat ) then

  echo "Did you remember to trim both ends of the file??"

  set origin=`head -1 trimmed.dat | awk '{print $2}'`

  awk '\!/#/{i+=((($2 - '$origin')/0.529177)*$3*2.5417462310548435); print $2, i}' trimmed.dat | tail -1

  endif

endif
