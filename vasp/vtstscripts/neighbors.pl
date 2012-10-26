eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

# 12-12-2002

#
# Calculates the distance, with periodic boundary condition,
# from a specific atom to all other atoms in the file.
#

#
# Open up and read the confile.
#

  $ift = index($ARGV[0],".con") + 1 ;
  unless ($ift) {die "First input argument must be a .con file" ;}

  open CON , $ARGV[0] ;
  while (<CON>) {$cf.=$_ ;}
  close CON ;
  @cf = split /\n/ , $cf ;

#
# Read box side lengths, the number of types of atoms
# and the number of each type.
#  
  
  $line = @cf[2] ; chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
  @box = @line[0..2] ;
  $line = @cf[6] ; chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
  $ntyp = @line[0] ;
  $line = @cf[7] ; chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
  @natoms = @line[0..$ntyp-1] ;

  $cumsum[0]=@natoms[0] ;
  for ($i=1; $i<$ntyp; $i++) {$cumsum[$i]+=$cumsum[$i-1] + $natoms[$i] ;}

#
# Calculate the distance from atom # $ARGV[1] to all other atoms and
# sort them in ascending order.
#

  $ca = @ARGV[1] ;
  @Rca = &coordinates($ca,@cf) ;

# Here we calculate the distances
  open OUT , ">out.tmp" ;
  $p = 0 ;
  for ($i=1; $i<=$cumsum[$ntyp-1]; $i++) {
    $k = $i + 10 + 2*$p ; 
    $line = @cf[$k] ; chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
    @Ra = @line[0..2] ; 
    @d = &dpc($Ra[0]-$Rca[0],$Ra[1]-$Rca[1],$Ra[2]-$Rca[2]) ;
    $dist = sqrt($d[0]**2 + $d[1]**2 + $d[2]**2) ;
    printf OUT "%4d %13.7f %13.7f %13.7f    ... %13.7f\n",$i,$Ra[0],$Ra[1],$Ra[2],$dist ;
#    print OUT $i,"   ",$Ra[0],"   ",$Ra[1],"   ",$Ra[2]," ...  ",$dist,"\n" ;
    if ($i == $cumsum[$p]) {$p++ ;}
  }
  close OUT ;

# Here we get the distances and index them.
  open IN , "out.tmp" ;
  $i=0 ;
  while ($line=<IN>) {
    $file.=$line ;
    chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
    $d[$i] = $line[5];
    $indx[$i] = $i ;
    $i++ ;
  }
  close IN ;
  unlink "out.tmp" ;

# Sort with selection sort.
  for ($i=0; $i<$cumsum[$typ-1]; $i++) {
    $iptr = $i ;
    for ($j=$i+1; $j<$cumsum[$typ-1]; $j++) {if($d[$j] < $d[$iptr]){$iptr = $j ;}}
    if ($i != $iptr) {
      $tmp = $d[$i] ;
      $d[$i] = $d[$iptr] ;
      $d[$iptr] = $tmp ;
      $tmp = $indx[$i] ;
      $indx[$i] = $indx[$iptr] ;
      $indx[$iptr] = $tmp ;
    }
  }

# Print out the sorted data.
  open OUT, ">neighdist.dat" ;
  @file = split /\n/ , $file ;
  for ($i=0; $i<$cumsum[$typ-1]; $i++) {
    print OUT $i,"   ",$file[$indx[$i]],"\n" ;}
  close OUT ;

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
    
  sub coordinates {
    my($in,@file) = @_ ;
    for ($i=0; $i<$ntyp; $i++) {
      if ($in <= $cumsum[$i]) {
        $line = @file[$in+10+2*$i] ; chomp($line) ; $line=~s/^\s+//g ;
        @line = split /\s+/,$line ;
        last ;
      }
    }
    @line[0..2] ;
  }

  sub dpc {
    my @d = @_ ;
    while($d[0] >  $box[0]/2)  {$d[0] -= $box[0] ;}
    while($d[0] < -$box[0]/2)  {$d[0] += $box[0] ;}
    while($d[1] >  $box[1]/2)  {$d[1] -= $box[1] ;}
    while($d[1] < -$box[1]/2)  {$d[1] += $box[1] ;}
    while($d[2] >  $box[2]/2)  {$d[2] -= $box[2] ;}
    while($d[2] < -$box[2]/2)  {$d[2] += $box[2] ;}
    @d ;
  }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
