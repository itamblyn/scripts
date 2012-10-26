#!/bin/tcsh

foreach machine ( selenium tellurium beryllium carbon neon nitrogen fluorine oxygen )

  echo -n $machine " "

  ssh ${machine} top -b -n1 | awk '/load/{print $10, $11, $12, "(5 min)", $14, "(15 min)"}' | sed s/,//g 

end
