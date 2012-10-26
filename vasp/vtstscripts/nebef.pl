eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

# Script prints the force, energy etc of OUTCAR's in immediate subdir
# of present working directory. Older versions of this script had specific
# paths hardcoded.  

  use Cwd ;
  $dir = cwd ;      
  @l1=`ls -la $dir/[0-9][0-9]/OUTCAR`;   #specifies location of outcars
  $i = 0 ;
  foreach $_ (@l1) {
    chop() ;
    @t1=split() ;
    $t2=$t1[@t1-1] ;
#    $steps  = `grep 'energy  without' $t2 | wc |cut -c 0-8` ; 
    $energy = `grep 'energy  without' $t2 | tail -n 1 |cut -c 68-78` ; 
    $force =  `grep 'max\ at' $t2 | tail -n 1 |cut -c 27-38` ;
#    chop($steps) ;
    chop($energy) ;
    chop($force) ;
    if(!$i) {$e0 = $energy ;}
    $rel = $energy - $e0 ;
    @f4 = ($i,$force,$energy,$rel) ;
    printf "%4i %16.8f %16.8f %16.8f \n",@f4 ; 

    $i++ ;

  };

