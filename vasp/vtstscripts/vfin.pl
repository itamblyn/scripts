eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";

# Check to see if there is an OUTCAR and read the ICHAIN variable from it
  if (-e "OUTCAR"){
    $_ = `head -1000 OUTCAR | grep -i ICHAIN ` ;
  }elsif (-e "01/OUTCAR"){
    $_ = `head -1000 01/OUTCAR | grep -i ICHAIN ` ;
  }else{
    die "No OUTCAR found !!! \n"
  }
  @s = split /\s+/ , $_ ;
  $ichain = $s[4] ;

  die "No result directory name given! \n" if @ARGV < 1 ;
  $dir = $ARGV[0] ;
  mkdir $dir, 0755 ;

# Clean up the run
  if ($ichain == 0){
# If there is a 00 directory then this is a NEB run
    if (-e "00") {
      system "$Bin/nebclean.sh $dir" ;
    }else{
      system "$Bin/vclean.sh $dir" ;
    }
  }elsif ($ichain == 1){
    system "$Bin/dymclean.sh $dir"
  }elsif ($ichain == 2){
    if((-e "MODECAR") && !(-e "01") && !(-e "02")){
      system "$Bin/dimclean.sh $dir" ;
    }elsif((-e "01") && (-e "02")){
      system "$Bin/dimclean2.sh $dir" ;
    }
  }elsif ($ichain == 3){
    system "$Bin/lanclean.sh $dir" ;
  }elsif ($ichain == 4){
    system "$Bin/insclean.sh $dir" ;
  }

