eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

$confile=@ARGV[0];
open IN, "$confile";
open OUT,">out.con";
$i=0;
while($line=<IN>){
  $i++;
  chomp($line);
  $line=~s/^\s+//g;
  if($i<=11) {print OUT $line,"\n";}
  @line=split /\s+/,$line;
  if($i==3){
    $ax=@line[0];
    $ay=@line[1];
    $az=@line[2];
    $ax_2=$ax/2;
    $ay_2=$ay/2;
    $az_2=$az/2;
  }
  if($i>11){    
    $line[0]-=$ax_2;
    @line[1]-=$ay_2;
    @line[2]-=$az_2;
    printf OUT "%13.10f %17.10f %17.10f %5d %5d \n",@line;
  }
}

