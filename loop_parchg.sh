#!/bin/tcsh


echo "VOLUME" > planechg.in
echo "3"    >> planechg.in
echo "0"    >> planechg.in
echo "0.0" >> planechg.in
echo "1" >> planechg.in
echo "0" >> planechg.in
echo " " >> planechg.in
echo "0" >> planechg.in
echo "0">> planechg.in


foreach r ( PARCHG* )

  rm VOLUME
  ln -s $r VOLUME
  set band=`echo $r | sed s/"\."/" "/g | awk '{print $2}'`
  ~/bin/planechg.x
  mv planechg.dat "$band".dat

end
