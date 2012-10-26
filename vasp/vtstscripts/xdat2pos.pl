eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

  die "PLEASE INPUT THE NUMBER OF THE STEP YOU WANT TO RETRIEVE!\n" if @ARGV < 1 ;

# Get the length of the XDATCAR file. This is bad but works for now !!!
  open TMP , "XDATCAR" ;
  $nl = 0 ;
  while(<TMP>){$nl++ ;} 
  close TMP ;

# Get information about, number and types of atoms, box lengths ets.
# from the POSCAR file and write out to the output POSCAR file
  open POS , "POSCAR" or die " NO POSCAR IN THIS DIRECTORY \n";
  open OUT , ">POSCAR.out" ; 
  for($i=0; $i<8; $i++){
    $pos = <POS> ; print OUT $pos ;
    chomp($pos) ; $pos=~s/^\s+//g ; @pos=split /\s+/,$pos ;
    if($i == 0){
      @elements = split /\s+/ , $pos ;
      $nel = @elements ;
    }
    if($i == 5){
      @not[0..$nel-1] = @pos[0..$nel-1] ;
      while($not[$k] != undef){$natoms+=$not[$k++] ;} 
    }
  }
  $pos = undef ;
  for($i=0; $i<$natoms; $i++){$pos .= <POS> ;}
  close POS ;
  @pos = split /\n/ , $pos ;

# Get the right step from the XDATCAR
  open XDAT , "XDATCAR" or die " NO XDATCAR IN THIS DIRECTORY \n";
  $a = $ARGV[0] ;
  $st = 6+($natoms+1)*($a-1) ;
  $fn = $st+$natoms ;  

  die "THE SEARCH IS OUT OF BOUNDS\n" if($st > $nl || $fn > $nl) ;

  for($i=0; $i<$st-1; $i++){$in = <XDAT> ;}
  $in = <XDAT> ;
  $j = 0 ; 
  for($i=$st+1; $i<=$fn; $i++){
    $in = <XDAT> ;
    if($in == undef){print "hallo\n" ;}
    chomp($in) ; $in=~s/^\s+//g ; @in=split /\s+/,$in ;
    $p = $pos[$j] ;
    $j++ ;
    chomp($p) ; $p=~s/^\s+//g ; @p=split /\s+/,$p ;
    printf OUT "%15.10f %15.10f %15.10f %5s %5s %5s \n",$in[0],$in[1],$in[2],$p[3],$p[4],$p[5] ;
#    $IN = <STDIN> ;
  }

  close OUT ;
  close XDAT ;

## Assign a type (element) to each atom.
#
#  $n = 0 ;
#  $j = 0 ;
#  for($i=1; $i<=$natoms; $i++){ 
#    $j++ ;
#    if($j <= $not[$n]){$type[$i-1] = $elements[$n] ;}
#    else{$n++ ; $type[$i-1] = $elements[$n] ; $j = 1 ;}  
##    print $i,"  ",$j,"  ",$type[$i-1];
##    $IN = <STDIN> ;
#  }
# 
## Read the XDAT file and make .xyz files
#
#  for($i=0; $i<5; $i++){$line = <XDAT> ; }                # Jump over the first few lines
#  $n=0 ;  
#  open MOV, ">movie.xyz" ;
#  while($line = <XDAT>){
#    chomp($line) ;
#    print MOV $natoms,"\n" ;
#    $f = $forces[$n] ; chomp($f) ; $f=~s/^\s+//g ; @f=split /\s+/,$f ;
#    $e = $energy[$n] ; chomp($e) ; $e=~s/^\s+//g ; @e=split /\s+/,$e ;
#    print MOV "FORCE:  $f[4]  ...  ENERGY:  $e[6]","\n" ;
#    for($i=0; $i<$natoms; $i++){
#      $line = <XDAT> ;
#      chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
##  Transform from direct coordinates to cart. coordinates.
#      $x = $line[0] ; $y = $line[1] ; $z = $line[2] ;
#      $xt=$x*$sidex[0]+$y*$sidey[0]+$z*$sidez[0];
#      $yt=$x*$sidex[1]+$y*$sidey[1]+$z*$sidez[1];
#      $zt=$x*$sidex[2]+$y*$sidey[2]+$z*$sidez[2];
#      $x=$latt*$xt; $y=$latt*$yt; $z=$latt*$zt;
#      printf MOV "%2s %18.13f %18.13f %18.13f \n",$type[$i],$x,$y,$z ;  
#    }
#    $n++ ;
#  }
