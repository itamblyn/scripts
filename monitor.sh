#!/bin/tcsh

foreach n (`wwnodes | awk '/itamblyn/{print $1}'`)

    echo -n $n " "

    ssh ${n} top -b -n1 | awk '/load/{print $10, $11, $12, $13, $14}' 

end
