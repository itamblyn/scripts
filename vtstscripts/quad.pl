eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

@args=@ARGV;
(@args==1||@args==2) || die "usage: quad.pl <inputfile> <outputfile>\n";
$inputfile=@args[0];
if(@args==2){$outputfile=@args[1];}
else{$outputfile="ciPOSCAR";}

open (IN,"<$inputfile");
open (OUT,">$outputfile");

for($line=0; $line<5; $line++){$out=<IN>; print OUT "$out";}

$natoms=<IN>;
$natoms=~s/^\s+//g;
@natoms=split(/\s+/,$natoms);
$totatoms=0;
for($i=0; $i<@natoms; $i++){
  $totatoms+=$natoms[$i];
  @natoms[$i]*=4;}
$natoms=join("  ",@natoms);
print OUT "$natoms\n";

for($line=6; $line<8; $line++){$out=<IN>; print OUT "$out";}

for($i=0; $i<$totatoms; $i++){
  $_=<IN>;
  $_=~s/^\s+//g;

  @out=split(/\s+/,$_);
  $out=join("  ",@out);
  print OUT "$out\n";
	
  @out=split(/\s+/,$_); 
  @out[0]-=1;
  $out=join("  ",@out);
  print OUT "$out\n";

  @out=split(/\s+/,$_); 
  @out[1]-=1;
  $out=join("  ",@out);
  print OUT "$out\n";

  @out=split(/\s+/,$_); 
  @out[0]-=1; @out[1]-=1;
  $out=join("  ",@out);
  print OUT "$out\n";
}

