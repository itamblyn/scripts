eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-
################################################################################
# get iteration information from OUTCAR file (or OUTCAR type file) 
################################################################################
  if (@ARGV==0) {
    $file='OUTCAR';
  } else {
    $file=$ARGV[0];
  };

  while (!-e $file) {
    print "\n\nfilename =? ";  chop($file = <STDIN>);
  };

  open(FILE, "grep Iteration $file |");
  @lines=<FILE>;
  close(FILE);

  $count=0;
  foreach $line (@lines) {
    $line=~m/(Iteration\s*([0-9]+)\s*\(\s*([0-9]+)\))/;
    $iterations[$count]=$1;
    $nuclear[$count]=$2;
    $electronic[$count]=$3;
    if ($count==$#lines) {
      print "$iterations[$count]\n";
    } elsif ($nuclear[$count] > $nuclear[$count-1]) {
      print "$iterations[$count-1]\n";
    };
    $count++;
  };
