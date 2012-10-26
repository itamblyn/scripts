eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

@args=@ARGV;
(@args==2) || die "usage: boxset.pl <POSCAR> <new_lattice_constant> \n";

$poscarfile1=$args[0];
($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description)
  =read_poscar($poscarfile1);

$karcoords=dirkar($coordinates,$basis,$lattice,$total_atoms);
$lattice_old=$lattice;
$lattice=$args[1];
for ($i=0;$i<3;$i++) {
  for ($j=0;$j<3;$j++) {
    $basis->[$i][$j]*=$lattice/$lattice_old;
  }
}
$coordsnew=kardir($karcoords,$basis,$lattice,$total_atoms);

print "Total atoms: $total_atoms...\n";
print "Lattice: $lattice...\n";
print "Shift: ".$shift[0]."  ".$shift[1]."  ".$shift[2]."\n";

write_poscar($coordsnew,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,"POSCAR.out");

