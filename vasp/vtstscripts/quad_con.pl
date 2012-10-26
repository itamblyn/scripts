eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

#-------------
# This script double the x and y box dimensions of the con file created by quad.pl.
# so that PBC can be applied to the created quadruple images.
#-------------
@args=@ARGV;
(@args==1||@args==2) || die "usage: quad_con.pl <inputfile> <outputfile>\n";
$inputfile=@args[0];
if(@args==2){$outputfile=@args[1];}
else{$outputfile="ciPOSCAR.con";}

open (IN,"<$inputfile");
open (OUT,">$outputfile");

for($line=0; $line<2; $line++){$out=<IN>; print OUT "$out";}

$box=<IN>;
$box=~s/^\s+//g;
@box=split(/\s+/,$box);
for($i=0; $i<2; $i++){
  @box[$i]*=2.0;}
$box=join("  ",@box);
print OUT "$box\n";

for($line=3; $line<6; $line++){$out=<IN>; print OUT "$out";}
$ntype=<IN>;
$ntype=~s/^\s+//g;
@ntype=split(/\s+/,$ntype);
$ntype=$ntype[0];
print OUT "$ntype\n";

$natoms=<IN>;
$natoms=~s/^\s+//g;
@natoms=split(/\s+/,$natoms);
$totatoms=0;
for($i=0; $i<@natoms; $i++){
  $totatoms+=$natoms[$i];
}
$natoms=join("  ",@natoms);
print OUT "$natoms\n";

$rest=8+$totatoms+2*$ntype+1;
for($line=8; $line<$rest; $line++){$out=<IN>; print OUT "$out";}
close(IN);
close(OUT);
