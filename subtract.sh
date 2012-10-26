#!/bin/csh


# subtract off energy of first configuration

set first=`head -1 energy.dat | awk '{print $2}'`
echo "# dist energy ab_energy" > energy_sub.dat
awk '{print $1, $2 - '$first', $2}' energy.dat >> energy_sub.dat

# calculate slope with numerical derivative

awk '{print $1, -($2 - n);n=$2}' energy.dat | tail +2 > difference.dat
