eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

@args=@ARGV;
@args==1 || die "usage: chgnorm.pl <reference CHGCAR> \n";


open (IN1,$args[0]) || die ("Can't open file $!");

for ($i=0;$i<5;$i++) {
   $line1=<IN1>;
   $header1.=$line1;
}

$atoms1=<IN1>;
$header1.=$atoms1;

@atoms1=split(/\s+/,$atoms1);

$sum1 += $_ for @atoms1;

print "Atoms in file1: ".$sum1."\n";

for ($i=0;$i<$sum1+2;$i++) {
   $header1.=<IN1>;
}

$points1=<IN1>;
$header1.=$points1;

@points1=split(/\s+/,$points1);

$psum1=1;

for ($i=1;$i<@points1;$i++) {
   $psum1*=$points1[$i];
}

$total = 0.0;
$counter = 0;

for ($i=0;$i<$psum1/5;$i++) {
   $line1=<IN1>;
   @line1=split(/\s+/,$line1);
   for ($j=1;$j<@line1;$j++) {
       $total = $total + $line1[$j];
       $counter = $counter + 1;
   }
}

print "\n TOTAL voxel value = ".($total)."\n";
print "\n AVERAGE voxel value = ".($total/$counter)."\n";
print "\n Number of voxels = ".($counter)."\n";

close(IN1);
