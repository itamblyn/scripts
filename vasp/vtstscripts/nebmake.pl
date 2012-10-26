eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

# 13-05-2005

# Makes a linear interpolation between to POSCARs with N points

use FindBin qw($Bin) ;
use lib "$Bin" ;
use Vasp ;

# Get the input parameters
  die "PLEASE INPUT THE TWO POSCARS AND THE NUMBER OF IMAGES!\n" if @ARGV < 3 ;
  $pos1 = $ARGV[0] ;
  $pos2 = $ARGV[1] ;
  $nim = $ARGV[2]+2 ;

  die "THE NUMBER OF IMAGES (INCLUDING END-POINTS) IS  LIMITED TO 100 AS IS \n" 
if $nim > 100 ;

# Read in the POSCAR files and make sure that the number of atoms and types are the 
# same
  ($coo1,$basis1,$lattice1,$natoms1,$totatoms1,$selectiveflag,$selective)
  =read_poscar($pos1) ;
  ($coo2,$basis2,$lattice2,$natoms2,$totatoms2,$selectiveflag,$selective)
  =read_poscar($pos2) ;
  if($totatom1 != $totatom2) {
    print "TOTAL NUMBER OF ATOMS IN FILE 1 IS:  ",$totatoms1,"\n" ;
    print "TOTAL NUMBER OF ATOMS IN FILE 2 IS:  ",$totatoms2,"\n" ;
    die ;
  }
  if(@{$natoms1} != @{$natoms2}) {
    print "TYPES OF ATOMS IN FILE 1 IS:  ",$n=@{$natoms1},"\n" ;
    print "TYPES OF ATOMS IN FILE 2 IS:  ",$n=@{$natoms2},"\n" ;
    die ;
  } 
  for($i=0 ; $i<@{$natoms1} ; $i++) {
    if($natoms1->[$i] != $natoms2->[$i]) {
      print "FOR ELEMENT ",$i," ... \n" ;
      print "... ATOMS IN FILE 1 IS:  ",$natoms1->[$i],"\n" ;
      print "... ATOMS IN FILE 2 IS:  ",$natoms2->[$i],"\n" ;
      die ;
    }
  } 
  if($lattice1 != $lattice2) {
    print "THE LATTICE CONSTANT IN FILE 1 IS:  ",$lattice1,"\n" ;
    print "THE LATTICE CONSTANT IN FILE 2 IS:  ",$lattice2,"\n" ;
    die ;
  }
  for($i=0 ; $i<3 ; $i++) {
    for($j=0 ; $j<3 ; $j++) {
      if($basis1->[$j][$i] != $basis1->[$j][$i]) {
        print "BASIS ELEMENT ",$i," ",$j," ... \n" ;
        print "... IS IN FILE 1:  ",$basis1->[$j][$i],"\n" ;
        print "... IS IN FILE 2:  ",$basis2->[$j][$i],"\n" ;
        die ;
      }
    }
  }

# Ok, the POSCARs appear to be for the same system.
# Get te header, i.e. the element symbols from the first POSCAR line 
  $header=`head -n 1 $pos1` ;
  chop($header) ;

# Calculate the distance between the two images ... dirkar: direct -> cartesian 
  $diff = pbc_difference($coo2,$coo1,$totatoms1) ;
  dirkar($diff,$basis1,$lattice1,$totatoms1) ;
  for($i=0 ; $i<$totatoms1 ; $i++) {
    for($j=0 ; $j<3 ; $j++) {
      $step->[$i][$j] = ($diff->[$i][$j])/($nim-1) ;
    }
  }
# Because zero is the number of the first image
  $nim-- ;

# Put the POSCAR in the initial state folder
  mkdir "00" ;
  write_poscar($coo1,$basis1,$lattice1,$natoms1,$totatoms1,
               $selectiveflag,$selective,$header,"00/POSCAR") ;  
# Put the POSCAR in the final state folder
  if($nim < 10) {$dir = "0$nim" ;}
  else {$dir = "$nim" ;}
  mkdir $dir ;
  write_poscar($coo2,$basis1,$lattice1,$natoms1,$totatoms1,
               $selectiveflag,$selective,$header,"$dir/POSCAR") ;

# Make the rest of the images in the chain
  for($i=0 ; $i<$totatoms1 ; $i++) {
    for($j=0 ; $j<3 ; $j++) {
      $t->[$i][$j] = $coo1->[$i][$j] ;
    }
  }
  for($im=1 ; $im<$nim ; $im++) {
    dirkar($t,$basis1,$lattice1,$totatoms1) ;
    for($i=0 ; $i<$totatoms1 ; $i++) {
      for($j=0 ; $j<3 ; $j++) {
#
# WHAT ABOUT ..PBC.. 
#
       $t->[$i][$j] += $step->[$i][$j] ;
      }
    }
    kardir($t,$basis1,$lattice1,$totatoms1) ;    
    if($im < 10) {$dir = "0$im" ;}
    else {$dir = "$im" ;}
    mkdir $dir ;
    write_poscar($t,$basis1,$lattice1,$natoms1,$totatoms1,
                 $selectiveflag,$selective,$header,"$dir/POSCAR") ;
  }

  print "\n" ;
  print "OK, ALL SETUP HERE\n" ;
  print "FOR LATER ANALYSIS, PUT OUTCARs IN FOLDERS 00 and " ;
  if($im < 10) {
    print "0$nim !!! \n" ;
  }else {
    print "$nim !!! \n" ;
  }

