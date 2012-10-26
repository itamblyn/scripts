eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

@ARGV>0 || die "usage: doslanalyze.pl <s,p,d,a=(all)> <atom num(s)>\n";

# If no orbital flag, plot the d-band
$colup=5;
$coldown=6;
if(@ARGV>0){
  $arg1=@ARGV[0];
  if($arg1 =~ /^\d+$/){
  }else{
    $oflag=shift(@ARGV);
    if($oflag eq 's'){$colup=1;$coldown=2;}
    if($oflag eq 'p'){$colup=3;$coldown=4;}
    if($oflag eq 'd'){$colup=5;$coldown=6;}
    if($oflag eq 'a'){$col=7} # analyze
}}

# Set selected atom(s)
$NumAtm=@ARGV;
@Atm=@ARGV;

# Find how many atom DOS there are
opendir MAINDIR, "." or die "can't open this dir!" ;
@DOSFILE = grep /^DOS\d+$/, readdir MAINDIR ;
$NumDOS=@DOSFILE;
closedir MAINDIR ;

# If no atom is selected, analyze all of them
if($NumAtm==0){ $NumDOS=$NumDOS-1; }
else{ $NumDOS=$NumAtm; }

$first=1;
# Read selected DOS files
for ($i=1; $i<=$NumDOS; $i++){
  if($NumAtm==0){ $DOSFILE="DOS"."$i"; }
  else{ $a=$Atm[$i-1];  $DOSFILE="DOS"."$a"; }

  $j=0;
  open (DOS,$DOSFILE);
  <DOS>;
  while ($in=<DOS>) {
    $in=~s/^\s+//g;
    @line=split(/\s+/,$in);
    $eneval=$line[0];
    if($colup<7){$dosval=$line[$colup]-$line[$coldown];}
    else{$dosval=$line[1]-$line[2]+$line[3]-$line[4]+$line[5]-$line[6];}
    if($first){
      $Ene[$j]=$eneval;
      $Dos[$j]=$dosval;
    }else{
      $Ene[$j]==$eneval || die "Energy $j in file $DOSFILE does not match first file.\n";
      $Dos[$j]+=$dosval;
    }
    $j++;
  }
  close(DOS);
  if($first){$numene=$j; $first=0;}
  else{$numene==$j || die "Number of energy values in $DOSFILE does not match first file.\n";}
}

# Calculate the average energy state
$dossum=$dstsum=$dsqsum=0;
for($i=0; $i<$numene; $i++){
  $dossum+=$Dos[$i];
  $dstsum+=$Ene[$i]*$Dos[$i];
  $dsqsum+=$Ene[$i]*$Ene[$i]*$Dos[$i]*$Dos[$i];
}
$dossum || die "Total DOS is zero.\n";
$eneavg=$dstsum/$dossum;
$enevar=($dstsum*$dstsum-$dsqsum)/($numene*$NumDOS*$dossum*$dossum);
$enestd=sqrt($enevar);

print "Total States: $dossum\n"; 
print "Average Energy: $eneavg\n";
print "Standard Deviation: $enestd\n";

