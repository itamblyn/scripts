eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

@args=@ARGV;
@args==2 || die "usage: chgavg.pl <reference CHGCAR> <CHGCAR2>\n";


open (IN1,$args[0]) || die ("Can't open file $!");
open (IN2,$args[1]) || die ("Can't open file $!");
open (OUT,">CHGCAR_avg");

for ($i=0;$i<5;$i++) {
   $line1=<IN1>;
   $line2=<IN2>;
   $header1.=$line1;
}

$atoms1=<IN1>;
$header1.=$atoms1;
$atoms2=<IN2>;

@atoms1=split(/\s+/,$atoms1);
@atoms2=split(/\s+/,$atoms2);

$sum1 += $_ for @atoms1;
$sum2 += $_ for @atoms2;

print "Atoms in file1: ".$sum1.", Atoms in file2: ".$sum2."\n";

for ($i=0;$i<$sum1+2;$i++) {
   $header1.=<IN1>;
}
for ($i=0;$i<$sum2+2;$i++) {
   $dummy=<IN2>;
}

$points1=<IN1>;
$header1.=$points1;
$points2=<IN2>;

@points1=split(/\s+/,$points1);
@points2=split(/\s+/,$points2);

$psum1=1;
$psum2=1;

for ($i=1;$i<@points1;$i++) {
   $psum1*=$points1[$i];
   $psum2*=$points2[$i];
}

print "Points in file1: ".$psum1.", Points in file2: ".$psum2."\n";

if ($psum1 != $psum2) {die ("Number of points not same in two files!");}

print OUT $header1;

for ($i=0;$i<$psum1/5;$i++) {
   $line1=<IN1>;
   $line2=<IN2>;
   @line1=split(/\s+/,$line1);
   @line2=split(/\s+/,$line2);
   for ($j=1;$j<@line1;$j++) {
       $line1[$j]=($line2[$j]+$line1[$j])/2;
   }
   printf OUT "%20.10e %20.10e %20.10e %20.10e %20.10e\n",$line1[1],$line1[2],$line1[3],$line1[4],$line1[5];  

}

close(OUT);
close(IN2);
close(IN1);
