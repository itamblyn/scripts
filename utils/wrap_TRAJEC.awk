awk '! /H/{print $0} /H/{L=15.031741;N=1000000;x=$2-(int($2/L+N+0.5)-N)*L;y=$3-(int($3/L+N+0.5)-N)*L;z=$4-(int($4/L+N+0.5)-N)*L;print "C",x,y,z;}' TRAJEC.xyz > TRAJEC_C_pbc.xyz
