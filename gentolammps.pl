#!/usr/bin/perl

open (DFTB, "<$ARGV[0]") || die "Can't open input: $ARGV[0]\n";
$ARGV[0] =~ s/.gen//;
$ARGV[0] = "data.".$ARGV[0];
open (LAMP, ">$ARGV[0]") || die "Can't open output: $ARGV[0]\n";
@line = split " ", <DFTB>;
$natom = $line[0];
print (LAMP "# Position data for H2O\n");
print (LAMP "\n");
print (LAMP "$natom atoms\n");
@types = split " ", <DFTB>;
$ntypes = @types;
print (LAMP "$ntypes atom types\n");
print (LAMP "\n");
for ($i = 0; $i < $natom; $i++) {
  @line = split " ", <DFTB>;
  $atom[$i] = $line[1];
  $xcoord[$i] = $line[2];
  $ycoord[$i] = $line[3];
  $zcoord[$i] = $line[4];
}
$line = <DFTB>;
@a1 = split " ", <DFTB>;
@a2 = split " ", <DFTB>;
@a3 = split " ", <DFTB>;
close (DFTB);
#$ham = -0.5*$a1[0];
#$hap =  0.5*$a1[0];
$ham = 0.0;
$hap =  $a1[0];
print (LAMP "$ham $hap xlo xhi\n");
#$ham = -0.5*$a2[1];
#$hap =  0.5*$a2[1];
$ham = 0.0;
$hap =  $a2[1];
print (LAMP "$ham $hap ylo yhi\n");
#$ham = -0.5*$a3[2];
#$hap =  0.5*$a3[2];
$ham = 0.0;
$hap = $a3[2];
print (LAMP "$ham $hap zlo zhi\n");
print (LAMP "\n");
print (LAMP "$a2[0] $a3[0] $a3[1] xy xz yz\n");
print (LAMP "\n");
print (LAMP "Masses\n");
print (LAMP "\n");
for ($i = 1; $i <= $ntypes; $i++) {
  if ($types[$i-1] eq "C") {
    $mass = 12.0000;
  } elsif ($types[$i-1] eq "H") {
    $mass = 1.0080;
  } elsif ($types[$i-1] eq "O") {
    $mass = 15.9990;
  } elsif ($types[$i-1] eq "N") {
    $mass = 14.0000;
  } elsif ($types[$i-1] eq "S") {
    $mass = 32.065;
  } elsif ($types[$i-1] eq "Al") {
    $mass = 27.0000;
  } elsif ($types[$i-1] =~ /Si/) {
    $mass = 28.0855;
  } elsif ($types[$i-1] eq "F") {
    $mass = 19.0000;
  } elsif ($types[$i-1] eq "CL") {
    $mass = 35.5000;
  } else {
    die "No mass for element: $types[$i-1]\n";
  }
  printf (LAMP "%d %lf\n", $i, $mass);
}
print (LAMP "\n");
print (LAMP "Atoms\n");
print (LAMP "\n");
for ($i = 1; $i <= $natom; $i++) {
  print (LAMP "     $i $atom[$i-1] 0     $xcoord[$i-1] $ycoord[$i-1] $zcoord[$i-1]\n");
}
print (LAMP "\n");
close (LAMP);
