eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

###########################################################################
#                                                                         # 
# Goes though all the image folders in a NEB                              #
# run and forms a xyz-file from the POSCAR.                               #
# Note that the first line in the POSCAR has                              #
# to have the elemental symbols in the same                               #
# order as in the POTCAR, i.e. for a run the                              #
# methyl and a Ni slab                                                    #
#                                                                         #
# Ni C H                                                                  #
#  4.97615083550000037                                                    # 
#      1.0000000000000000    0.0000000000000000    0.0000000000000000     #
#      0.0000000000000000    0.8660254037832008    0.0000000000000000     #
#      0.0000000000000000    0.0000000000000000    3.6172536956853296     #
#   16   1   4                                                            #
# Selective dynamics                                                      #
# Direct                                                                  # 
#                                                                         #
###########################################################################

  use Cwd;
  use FindBin qw($Bin);
  use lib "$Bin";

  if(!$ARGV[0]){
    print "\nUsing POSCARs to generate movie\n\n" ;
    $filetype = "POSCAR" ;}
  else{
    print "\nUsing CONTCARs to generate movie\n\n" ;
    $filetype = "CONTCAR" ;
  }

  $dir = cwd;
  $zip = $ENV{'VTST_ZIP'} ;
  if($zip eq ''){ $zip = 'gzip' ; }
  $i=0;
  $string = "00" ;
  while(chdir $string){
# Grab the forces and energies to add to the xyz files
    $zipped = 0 ;
    $outcar = 1 ;
    if(-e "OUTCAR") { ; }
    elsif(-e "OUTCAR.gz"){ $zipped = 1 ; system "gunzip OUTCAR.gz" ; }
    elsif(-e "OUTCAR.bz2"){ $zipped = 1 ; system "bunzip2 OUTCAR.bz2" ; }
    else{ $outcar = 0 ; }

    if($outcar){
      $for = `grep 'Forces: m' OUTCAR | tail -1` ;
      if($for == undef){$for = `grep 'FORCES: m' OUTCAR | tail -1` ;}
      $ene = `grep 'energy  w' OUTCAR | tail -1` ;
      $f = $for ; chomp($f) ; $f=~s/^\s+//g ; @f=split /\s+/,$f ;
      $e = $ene ; chomp($e) ; $e=~s/^\s+//g ; @e=split /\s+/,$e ;
      if($i == 0){$e0 = $e[6] ;}
      if($zipped){system "$zip -9 OUTCAR &" ;}
    }

    $if = 1 ;
    if(!(-e $filetype)){print "copying\n"; system "cp POSCAR CONTCAR" ; $if = 0 ;}
    system "vasp2con.pl $filetype POSCAR.con > /dev/null" ;
    system "con2xyz.pl POSCAR.con > /dev/null" ;
    unlink "POSCAR.con" ;
    if(!$if){print "unlinking\n"; unlink "CONTCAR"}
    if($outcar){
      $e = $e[6] - $e0 ;
      system "sed s/'Generated with con2xyz'/'F:  $f[4]  ...  E:  $e'/g POSCAR.xyz > POSCAR_TMP.xyz" ;
    }else{
      system "sed s/'Generated with con2xyz'/' NO OUTCAR FOUND '/g POSCAR.xyz > POSCAR_TMP.xyz" ;
    }
    system "cp POSCAR_TMP.xyz ..//p$i.xyz" ;
    unlink "POSCAR.xyz" ;
    rename "POSCAR_TMP.xyz" , "POSCAR.xyz" ;

    $i++;
    if($i < 10) {$string = "0$i" ;}
    elsif($i < 100) {$string = "$i" ;}
    else {die "Too many images" ;}
    chdir $dir 
  }
  
  if (-e movie.xyz){
    unlink "movie.xyz" ;} 
  $i--;

  if($i < 10){
    system "cat p[0-9].xyz > movie.xyz" ;}
  elsif($i < 100){
    system "cat p[0-9].xyz p[0-9][0-9].xyz > movie.xyz" ;}
  else{print " TOO MANY XYZ FILES ...\n" ;}

#  system "rm -f [0-9][0-9]/POSCAR.xyz" ;

  unlink glob "p*.xyz"

