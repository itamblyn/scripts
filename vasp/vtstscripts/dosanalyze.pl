eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

@ARGV==0 || die "usage: dosanalyze.pl\n";

$j=0;
open (DOS,"DOS0");
while ($in=<DOS>) {
  $in=~s/^\s+//g;
  @line=split(/\s+/,$in);
  $Ene[$j]=$line[0];
  $ene=$line[0];
  $Dos[$j]=$line[1];
  $dos=$line[1];
  $j++;
}
$numene=$j;
close(DOS);

# Calculate the average energy state
$dossum=$dstsum=$dsqsum=0;
for($i=0; $i<$numene; $i++){
  $dossum+=$Dos[$i];
  $dstsum+=$Ene[$i]*$Dos[$i];
  $dsqsum+=$Ene[$i]*$Ene[$i]*$Dos[$i]*$Dos[$i];
}
$dossum || die "Total DOS is zero.\n";
$eneavg=$dstsum/$dossum;
$enevar=($dstsum*$dstsum-$dsqsum)/($numene*$dossum*$dossum);
$enestd=sqrt($enevar);

print "Total States: $dossum\n"; 
print "Average Energy: $eneavg\n";
print "Standard Deviation: $enestd\n";

