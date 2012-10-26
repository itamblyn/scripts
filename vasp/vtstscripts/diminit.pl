eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

#----------------------------------------------------------------------
# This script gets a POSCAR and a DISPLACECAR_sp to generate single-image 
# dimers for vasp dimer calculations.
#
# Last Modified Jan. 29, 2007 by Lijun Xu and Graeme Henkelman, UT-Austin
#----------------------------------------------------------------------

use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

(@ARGV>0) || die " USAGE: diminit.pl DIRECTORY OR NUMBER (OPTIONAL: DisplaceAlgo DisplaceRange Rcut MaxCord POSCAR DISPLACECAR_sp)\n"; 

$ndimers=0;
$DisplaceAlgo=0; #0 any atom allowed 1 lowest coordianted 2 less than a specified value
$DisplaceRange=0.01; # default: displace itself only
$rcut=0.01; # default: no neighbors, just itself
$MaxCord=8;
$poscarfilename="POSCAR";
$displacecarfilename="DISPLACECAR_sp";
# The distance between dimer images is now controlled by INCAR
# $deltaR=0.005 ; #0.01 is the distance between two images

if($ARGV[0]=~/^\d+$/){
  # It's a number, so we create that many dimers in the folders: pr000[1-n]
  $ndimers=$ARGV[0]; }
if(@ARGV>1){ $DisplaceAlgo=$ARGV[1]; }
if(@ARGV>2){ $DisplaceRange=$ARGV[2]; }
if(@ARGV>3){ $rcut=$ARGV[3]; }
if(@ARGV>4){ $MaxCord=$ARGV[4]; }
if(@ARGV>5){ $poscarfilename=$ARGV[5]; }
if(@ARGV>6){ $displacecarfilename=$ARGV[6]; }
#----------------------------------------------------------------------
#  Read in DISPLACEMENT file to decide which atom should be displaced
#  DISPALCEMENT has the same file structure as DISPLACECAR
#  Example:
#    0   0   0   1  : atom 1 will not be displaced
#    0.1 0.1 0.1 2  : atom 2, x, y, z will be displaced by a
#                     0.1-width Gaussian random number
#----------------------------------------------------------------------

open(DISPLACEMENT, "<$displacecarfilename")
  or die (" Error: cannot open the displacement file: ".$displacecarfilename."\n Note: default is DISPLACECAR_sp\n");
close(DISPLACEMENT);
($sigma, $total_atoms_disp) = read_othercar($displacecarfilename);
$DisplaceList=();
for ($i=0; $i<$total_atoms_disp; $i++) {
  #print $sigma->[$i][0]."  ".$sigma->[$i][1]."  ".$sigma->[$i][2]."  ".($i+1)."\n"; 
  if($sigma->[$i][0] != 0 || $sigma->[$i][1] != 0 || $sigma->[$i][2] != 0) {
    push @$DisplaceList, $i;
  }
}
if(@$DisplaceList==0){die "Error from diminit.pl; no atoms will be displaced. what is going on here. Check it.\n";}
#----------------------------------------------------------------------
#  Read POSCAR files and prepare  the random displacement 
#----------------------------------------------------------------------
print " Reading POSCAR\n";
($coordinates,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description)=read_poscar($poscarfilename);
($total_atoms_disp == $total_atoms) ||
   die " Error: Displacement file $displacecarfilename should have the same number of lines as there are atoms in the POSCAR file.\n";
set_bc($coordinates,$total_atoms);
$coordinates = dirkar($coordinates,$basis,$lattice,$total_atoms);
# build a list of atoms that be a center of displacement(default: displace all allowed to move) and their neighbors
($NcList,$NN_DisplaceList)=BuildDisplaceList($coordinates,$basis,$DisplaceList,$DisplaceAlgo,$DisplaceRange,$rcut,$MaxCord);
print "**********************************\n";
print "@$NcList\n";
print "**********************************\n";
for $i (@$DisplaceList){
  print "$i: @{$NN_DisplaceList->{$i}}\n";
}
#----------------------------------------------------------------------
#  Do the displacement to generate dimer(s)
#----------------------------------------------------------------------
if($ndimers==0){
  $NumSearches=1;
  $folder=$ARGV[0];
}else{
  $NumSearches=$ndimers;
}
print " Generating files for $NumSearches dimer searches\n";
srand();
for($k=1; $k<=$NumSearches; $k++){
  $randomnum=rand();
  print "randomnum=$randomnum\n";
  $atoms2bdisplaced=BuildNewDisplacecar($DisplaceList,$NcList,$NN_DisplaceList,$DisplaceAlgo,$randomnum);
  # folder name
  if($ndimers>0){
    $prnum=sprintf "%04d",$k ;
    $folder="pr".$prnum;
  }else{
    $folder=$ARGV[0]; }
  print " Dir: $folder\n";

  # position of dimer
  for($i=0;$i<$total_atoms;$i++) {
    for ($j=0; $j<3; $j++) {
      $coordinatesnew->[$i][$j]=$coordinates->[$i][$j];
    }
  }
  for $i (@$atoms2bdisplaced) {
    for ($j=0; $j<3; $j++) {
      $coordinatesnew->[$i][$j]=$coordinatesnew->[$i][$j]+$sigma->[$i][$j]*gauss;
    }
  }
  #write_othercar($coordinatesnew,$total_atoms,"random");

  # spring-relaxation: prevent atoms from being too close 
  print "relaxing the initial guess ... ... \n";
  spring_relaxation($coordinatesnew,$basis,$lattice,$total_atoms,$selective);
  print "... done\n";

  # displacement direction and magnitude
  $coordinatesnew = kardir($coordinatesnew,$basis,$lattice,$total_atoms);
  set_bc($coordinatesnew,$total_atoms);
  $coordinates = kardir($coordinates,$basis,$lattice,$total_atoms);
  $r_vector = pbc_difference($coordinatesnew,$coordinates,$total_atoms);
  $r_vector = dirkar($r_vector,$basis,$lattice,$total_atoms);
  $r_mag = magnitude($r_vector,$total_atoms);
  ($r_mag>0) || die " Error in diminit.pl: No atoms have been displaced.\n";
  $reciprocal= 1.0/$r_mag;
  $r_vector = vmult($r_vector,$reciprocal,$total_atoms);

  # write dimer images
  if(-e $folder) { die " Error from diminit.pl: $folder already exists.\n"; }
  system "mkdir $folder";
  $outputfilename="ciPOSCAR";
  write_poscar($coordinatesnew,$basis,$lattice,$num_atoms,$total_atoms,$selectiveflag,$selective,$description,$outputfilename);
  system "mv ciPOSCAR $folder/POSCAR";
  write_othercar($r_vector,$total_atoms,$outputfilename);
  system "mv ciPOSCAR $folder/MODECAR";
  system "cp KPOINTS POTCAR $folder/";
  $coordinates = dirkar($coordinates,$basis,$lattice,$total_atoms);

  # fix this later.
  if(-e "INCAR") { system "cp INCAR $folder/"; }
  if(-e "akmc.sub") { system "cp akmc.sub $folder/"; }
}

# ---------------------------------------------------------------------------------------------------------
# get all the candidate atoms and the neighborlist between them: ($NcList,$NN_DisplaceList)=
# BuildDisplaceList($R,$basis,$DisplaceList,$DisplaceAlgo,$DisplaceRange,$rcut,$MaxCord);
# ---------------------------------------------------------------------------------------------------------
sub BuildDisplaceList{
  my ($R,$basis,$DisplaceList,$DisplaceAlgo,$DisplaceRange,$rcut,$MaxCord)=@_;
  my ($vector1,$vector2,$vector3,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_);
  my ($NN_NcordFile,$NN,$NCord,$i,$j,$good,$NcList,$NN_DisplaceList,$dummy);
  $NcList={};
  $NN_DisplaceList={};
  $NN_NcordFile="neighborlist.dat";
  $good=0;
  $j=@$R;
  print "Rcut=$rcut\n";
  if(-e $NN_NcordFile){
    ($NN,$NCord,$good)=readNN_Ncord($NN_NcordFile,$rcut,$j);
  }
  if(!$good){ # need to make a new list
    for($j=0;$j<3;$j++) {
      $vector1->[0][$j]=$basis->[$j][0];
      $vector2->[0][$j]=$basis->[$j][1];
      $vector3->[0][$j]=$basis->[$j][2];
    }     
    $Ax=magnitude($vector1,1);
    $Ay=magnitude($vector2,1);
    $Az=magnitude($vector3,1);
    print "Ax=$Ax Ay=$Ay Az=$Az\n";
    $Ax2=0.5*$Ax;
    $Ay2=0.5*$Ay;
    $Az2=0.5*$Az;
    $Ax2_=-0.5*$Ax;
    $Ay2_=-0.5*$Ay;
    $Az2_=-0.5*$Az;
    ($NN,$NCord)=FindNeighbors($R,$rcut,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_);
    print "**********************************\n";
    for $j (sort {$a <=> $b} keys %$NN){
      print $j."  ".$NCord->{$j}."  @{$NN->{$j}}\n";
    }
    writeNN_Ncord($NN_NcordFile,$NN,$NCord,$rcut);
  }
  $NN_DisplaceList=FindDLNeighbors($R,$DisplaceList,$DisplaceRange,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_);
  SWITCH: {   
   if($DisplaceAlgo==0){$NcList=FindLowestCoordAtom($DisplaceList,$NCord); last SWITCH;}
   if($DisplaceAlgo==1){$NcList=FindLessMaxcord($DisplaceList,$NCord,$MaxCord);last SWITCH;}
   if($DisplaceAlgo==2 || $DisplaceAlgo==3){
    $NcList=FindAtomsInIslands($R,$DisplaceList,$NN,$NCord,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_,$rcut,$DisplaceAlgo,$MaxCord);
    last SWITCH;
   }
   $NcList=$DisplaceList; #default: select a random atom
  }
  return ($NcList,$NN_DisplaceList);
}

sub FindNeighbors{
  my ($R,$rcut,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_)=@_;
  my ($Rijx,$Rijy,$Rijz,$Rij,$i,$j,$total_atoms);
  my ($NCord,$NN,@dummy);
  $total_atoms=@$R-1; 
  print "rcut=$rcut\n";
  for($i=0;$i<=$total_atoms;$i++){
    $NCord->{$i}=0;
    $NN->{$i}=[()];
  }
  for($i=0;$i<$total_atoms;$i++){
    for($j=$i+1;$j<=$total_atoms;$j++){
      $Rijx=$R->[$i][0] - $R->[$j][0];
      $Rijy=$R->[$i][1] - $R->[$j][1];
      $Rijz=$R->[$i][2] - $R->[$j][2];
      if($Rijx < $Ax2_){$Rijx+=$Ax;}
      if($Rijy < $Ay2_){$Rijy+=$Ay;}
      if($Rijz < $Az2_){$Rijz+=$Az;}
      if($Rijx > $Ax2){$Rijx-=$Ax;}
      if($Rijy > $Ay2){$Rijy-=$Ay;}
      if($Rijz > $Az2){$Rijz-=$Az;}
      $Rij=sqrt($Rijx*$Rijx+$Rijy*$Rijy+$Rijz*$Rijz);
      if($Rij < $rcut){
        $NCord->{$i}++;
        push @{$NN->{$i}}, $j;
        $NCord->{$j}++;
        push @{$NN->{$j}}, $i;
      }
    }
  }
  @dummy=sort {$a <=> $b} values %$NCord;
  $NCord->{"mincord"}=$dummy[0];
  return ($NN,$NCord);
}

sub FindDLNeighbors{
  my ($R,$DisplaceList,$rcut,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_)=@_;
  my ($Rijx,$Rijy,$Rijz,$Rij,$i,$j,$m,$n,$total_atoms,$NN);
  $total_atoms=@$DisplaceList-1;
  print "rcut=$rcut\n";
  $NN={};
  for $i (@$DisplaceList){
    $NN->{$i}=[()];
  }
  for($m=0;$m<$total_atoms;$m++){
    $i=$DisplaceList->[$m];
    for($n=$m+1;$n<=$total_atoms;$n++){
      $j=$DisplaceList->[$n];
      $Rijx=$R->[$i][0] - $R->[$j][0];
      $Rijy=$R->[$i][1] - $R->[$j][1];
      $Rijz=$R->[$i][2] - $R->[$j][2];
      if($Rijx < $Ax2_){$Rijx+=$Ax;}
      if($Rijy < $Ay2_){$Rijy+=$Ay;}
      if($Rijz < $Az2_){$Rijz+=$Az;}
      if($Rijx > $Ax2){$Rijx-=$Ax;}
      if($Rijy > $Ay2){$Rijy-=$Ay;}
      if($Rijz > $Az2){$Rijz-=$Az;}
      $Rij=sqrt($Rijx*$Rijx+$Rijy*$Rijy+$Rijz*$Rijz);
      if($Rij < $rcut){
        push @{$NN->{$i}}, $j;
        push @{$NN->{$j}}, $i;
      }
    }
  }
  return $NN;
}

sub writeNN_Ncord{
  my ($NN_NcordFile,$NN,$NCord,$rcut)=@_;
  my ($i,$mincord_total);
  open(ISMIN,">$NN_NcordFile") || die "cannot open $NN_NcordFile\n";
  print ISMIN "Rcut  ".$rcut."\n";
  print ISMIN "mincord  ".$NCord->{"mincord"}."\n";
  print ISMIN "atom\t"."ncoord\t"."neighbors\n";
  for $i (sort {$a <=> $b} keys %$NN){
    print ISMIN $i."\t".$NCord->{$i}."\t"."@{$NN->{$i}}\n";
  }
  close ISMIN;
  if(($i+1)==$total_atoms){
    die "Error in writeNN_Ncord: # of neighborlist elements is not same as # of atoms from POSCAR\n";
  }
}

sub readNN_Ncord{
  my ($NN_NcordFile,$rcut,$total_atoms)=@_;
  my ($NCord,$NN,@line,$line,$i,$j,$numatoms,$good);
  $NN={};
  $NCord={};
  $i=$j=$good=0;
  open(ISMIN,"<$NN_NcordFile") || die "cannot open $NN_NcordFile\n";
  while($line=<ISMIN>){
    $line=~s/^\s+//;
    if($line eq "\n"){next;}
    @line=split(/\s+/,$line);
    if(lc($line[0]) eq "rcut"){
      if($line[1]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/ && $line[1]==$rcut){
        $good =1;
        $i++;
        next;
      }else{
        last;
      }
    }
    if($i==1 && lc($line[0]) eq "mincord"){
      if($line[1]=~/^\d+$/) { $NCord->{"mincord"}=$line[1]; }
      $i++;
      next;
    }
    if($i==2 && lc($line[0]) eq "atom"){
      $i++;
      next;
    }
    if($i==3){
      $numatoms=@line-1;
      $NCord->{$line[0]}=$line[1];
      $NN->{$line[0]}=[(@line[2 .. $numatoms])];
      $j++;
    }
  }
  close ISMIN;
  if($j && $j != $total_atoms){$good=0;}
  return ($NN,$NCord,$good);
}

sub FindLowestCoordAtom{
  my ($DisplaceList,$NCord)=@_;
  my $NcList=();
  my ($Nc,$i,$mincord);
  $mincord=$NCord->{$DisplaceList->[0]};
  for $i (@$DisplaceList) {
    $Nc=$NCord->{$i};
    if($Nc < $mincord){
      $mincord=$Nc;
      $NcList=();
      push @$NcList, $i;
    }elsif($Nc == $mincord){
      push @$NcList, $i;
    }
  }
  if($NCord->{"mincord"} < $mincord){
    print "warning; the overall lowest coordinated atom(s) not included in the displacement. Make sure that is what you want.\n";
  }
  return $NcList;
}

sub FindLessMaxcord{
  my ($DisplaceList,$NCord,$MaxCord)=@_;
  my $NcList=();
  my $i=0;
  for $i (@$DisplaceList) {
    if($NCord->{$i} < $MaxCord){
      push @$NcList, $i;
    }
  }
  if(@$NcList == 0){
    print "MaxCord is too small($MaxCord) and no atoms in DISPLACECAR_sp are selected. Do random selection\n";
    $NcList=$DisplaceList;
  }
  return $NcList;
}

# ---------------------------------------------------------------------------------------------------------
# get a list of atoms to be displaced :  $atoms2bdisplaced=
# BuildNewDisplacecar($DisplaceList,$NcList,$total_atoms,$DisplaceAlgo,$DisplaceRange,$rcut);
# ---------------------------------------------------------------------------------------------------------
sub BuildNewDisplacecar{
  my ($DisplaceList,$NcList,$NN_DisplaceList,$DisplaceAlgo,$randomnum)=@_;
  my ($atoms2bdisplaced,$Nc,$atom);
  #pick out one displacing mechanism, currently it looks trivial
  SWITCH: {   
   if($DisplaceAlgo==0){$Nc=SelectNcAtom($NcList,$DisplaceAlgo,$randomnum);last SWITCH;}
   if($DisplaceAlgo==1){$Nc=SelectNcAtom($NcList,$DisplaceAlgo,$randomnum);last SWITCH;}
   if($DisplaceAlgo==2){$Nc=SelectNcAtom($NcList,$DisplaceAlgo,$randomnum);last SWITCH;}
   $Nc=SelectNcAtom($NcList,$DisplaceAlgo,$randomnum);  # default: displace around a random atom in DisplaceList
  }
  $atoms2bdisplaced=[($Nc)];
  for $atom (@{$NN_DisplaceList->{$Nc}}){ push @$atoms2bdisplaced, $atom;}
  print "Nc=$Nc leads to atoms to be displaced:  "."@$atoms2bdisplaced\n";
  return  $atoms2bdisplaced;
}

sub SelectNcAtom{
  my ($NcList,$DisplaceAlgo,$randomnum)=@_;
  my ($i,$j,$Nc,$total);
  print "randomnum=$randomnum\n"; 
  $j=@$NcList;
  $j=int($j*$randomnum);
  print "j=$j\n";
  for($i=0;$i<@$NcList;$i++){
    if(!$j){
      $Nc=$NcList->[$i];last;
    }else{ $j--;}
  }
  return $Nc;
}

# ---------------------------------------------------------------------------------------------------------
# relax an initial saddle point guess to prevent atoms from getting too close to each other
# spring_relaxation($coordinates,$basis,$lattice,$totalatoms,$selective);
# ---------------------------------------------------------------------------------------------------------
sub spring_relaxation{
  my $R=shift;
  my $basis=shift;
  my $lattice=shift;
  my $totalatoms=shift;
  my $selective=shift;
  my $epsilon=0.1;
  my $cutoff=0.3;
  my $stepmax=100000;
  my $drmax=0.2;
  my @line=();
  my ($i,$j,$Rijx,$Rijy,$Rijz,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_,$Rij);
  my ($force,$vector1,$vector2,$vector3,$fij,$phiij,$step,$goodsound,$converged,$frozen,$dr);
  # figure out box dimension [note: $R must be Cartesian]
  for($j=0;$j<3;$j++) {
    $vector1->[0][$j]=$basis->[$j][0];
    $vector2->[0][$j]=$basis->[$j][1];
    $vector3->[0][$j]=$basis->[$j][2];
  }
  $Ax=magnitude($vector1,1);
  $Ay=magnitude($vector2,1);
  $Az=magnitude($vector3,1);
  $Ax2=0.5*$Ax;
  $Ay2=0.5*$Ay;
  $Az2=0.5*$Az;
  $Ax2_=-0.5*$Ax;
  $Ay2_=-0.5*$Ay;
  $Az2_=-0.5*$Az;
  $Ax_=-1.0*$Ax;
  $Ay_=-1.0*$Ay;
  $Az_=-1.0*$Az;
  # figure out frozen flags for each degree of freedom
  for($i=0;$i<$totalatoms;$i++) {
    @line=split(/\s+/,$selective->[$i]);
    for($j=0;$j<3;$j++) {
      $frozen->[$i][$j]=$line[$j];
    }
  }
  # start the steepest descent loop
  $converged=0;
  for($step=0;$step<$stepmax;$step++) {
    if($converged) {last;}
    #print "Step No. $step\n";
    for($i=0;$i<$totalatoms;$i++) {
      for($j=0;$j<3;$j++) {
        $force->[$i][$j]=0.0;
      }
    }
    $goodsound=1;
    for($i=0;$i<$totalatoms-1;$i++) {
      for($j=$i+1;$j<$totalatoms;$j++) {
        $Rijx=$R->[$i][0] - $R->[$j][0];
        $Rijy=$R->[$i][1] - $R->[$j][1];
        $Rijz=$R->[$i][2] - $R->[$j][2];
        if($Rijx < $Ax2_){$Rijx+=$Ax;}
        if($Rijy < $Ay2_){$Rijy+=$Ay;}
        if($Rijz < $Az2_){$Rijz+=$Az;}
        if($Rijx > $Ax2){$Rijx-=$Ax;}
        if($Rijy > $Ay2){$Rijy-=$Ay;}
        if($Rijz > $Az2){$Rijz-=$Az;}
        $Rij=sqrt($Rijx*$Rijx+$Rijy*$Rijy+$Rijz*$Rijz);
        if($Rij < $cutoff){
          $goodsound=0;
          ($fij,$phiij)=linear_repulsion_pot($Rij,$cutoff);
          #($fij,$phiij)=exp_repulsion_pot($Rij,$cutoff);
          $force->[$i][0]=$force->[$i][0]-$fij*$Rijx/$Rij;
          $force->[$i][1]=$force->[$i][1]-$fij*$Rijy/$Rij;
          $force->[$i][2]=$force->[$i][2]-$fij*$Rijz/$Rij;
          $force->[$j][0]=$force->[$j][0]+$fij*$Rijx/$Rij;
          $force->[$j][1]=$force->[$j][1]+$fij*$Rijy/$Rij;
          $force->[$j][2]=$force->[$j][2]+$fij*$Rijz/$Rij;
        }
      }
    }
    if($goodsound) {
      $converged=1;
      next;
    }
    for($i=0;$i<$totalatoms;$i++) {
      for($j=0;$j<3;$j++) {
        #if($frozen->[$i][$j] == 1) {next;}
        if(lc($frozen->[$i][$j]) eq "f") {next;}
        $dr=$epsilon*$force->[$i][$j];
        if(abs($dr) > $drmax) {$dr=$drmax*$dr/abs($dr);}
        #print "$i: $j -- dr: $dr\n"; 
        $R->[$i][$j]+=$dr;
      }
      while ($R->[$i][0] < 0) {$R->[$i][0]+=$Ax;}
      while ($R->[$i][1] < 0) {$R->[$i][1]+=$Ay;}
      while ($R->[$i][2] < 0) {$R->[$i][2]+=$Az;}
      while ($R->[$i][0] > $Ax) {$R->[$i][0]-=$Ax;}
      while ($R->[$i][1] > $Ay) {$R->[$i][1]-=$Ay;}
      while ($R->[$i][2] > $Az) {$R->[$i][2]-=$Az;}
    }
  }
  if($converged == 1) {
    print "fully relaxed after $step (1: trivial) steps\n";
  }else {
    print "Not fully relaxed after $stepmax steps, but we can;t wait. move on\n";
  }
}

sub linear_repulsion_pot{ # E=kx-b, k< 0
  my $r=shift;
  my $cutoff=shift;
  my $slope=-0.3;
  my ($energy, $fij);
  my $intercept=$slope*$cutoff;
  $energy=$slope*$r-$intercept;
  $fij=$slope;
  return ($fij,$energy);
} 

sub exp_repulsion_pot{ # E=exp(kx)-b,k<0
  my $r=shift;
  my $cutoff=shift;
  my $slope=-0.3;
  my ($energy, $fij);
  my $intercept=exp($cutoff*$slope);
  my $item=exp($slope*$r);
  $energy=$item-$intercept;
  $fij=$slope*$item;
  return ($fij,$energy);
} 

# ---------------------------------------------------------------------------------------------------------
# divide adsorbates on the surface into islands and choose atoms from each island: $NcList=FindAtomsInIslands
# ($R,$DisplaceList,$NN,$NCord,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_,$rcut,$DisplaceAlgo,$MaxCord)
# ---------------------------------------------------------------------------------------------------------
sub FindAtomsInIslands{
  my ($R,$DisplaceList,$NN,$NCord,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_,$rcut,$DisplaceAlgo,$MaxCord)=@_;
  my ($islandsinfo,$islands,$NumIslands,$i,$j,$Nc,$total_atoms,$good,$mincord_island,$NcList);
  $islandsinfo="islands.dat";
  $total_atoms=@$DisplaceList;
  if(-e $islandsinfo){
    ($islands,$mincord_island,$good)=readIslands($islandsinfo,$rcut,$total_atoms);
  }
  if(!$good){ # create a new island information
    print "create a new island/coordination list\n";
    $islands=FindSurfIslands($R,$DisplaceList,$rcut,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_);
    $mincord_island=();
    $NumIslands=@$islands;
    for($i=0;$i<$NumIslands;$i++){      
      $mincord_island->[$i]=FindMostUndercoordianted($i,$islands,$NCord);
      print "isalnd $i: @{$islands->[$i]}\n";
      print "mincord in island $i is ".$mincord_island->[$i]."\n";
    }
    for($i=0;$i<@$DisplaceList;$i++){
      $j=$DisplaceList->[$i];
      print "atom $j: NCord=".$NCord->{$j}." @{$NN->{$j}}\n";
    }
    writeIslands($islandsinfo,$islands,$mincord_island,$rcut);
  }else{
    print "read in old islands and minimumcord information from file $island_minncord\n";
  }
  $NcList=();
  $NumIslands=@$islands;
  for($i=0;$i<$NumIslands;$i++){
    for $j (@{$islands->[$i]}){
      if($DisplaceAlgo==2){
        if($NCord->{$j}==$mincord_island->[$i]){ push @$NcList,$j;}
      }elsif($DisplaceAlgo==3){
        if($NCord->{$j} < $MaxCord){ push @$NcList,$j;}
      }
    }
  }
  if(@$NcList==0){
    print "MaxCord is too small($MaxCord) in FindAtomsInIslands. Have to choose atoms randomly\n";
    $NcList=$DisplaceList;}
  return $NcList;
}

sub FindSurfIslands{
  my ($R,$DisplaceList,$rcut,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_)=@_;
  my ($Rijx,$Rijy,$Rijz,$Rij);
  my ($i,$j,$m,$domains,$mark,$islands,$natoms);
  $domains={};
  $mark=0;
  for $i (@$DisplaceList){
    if(!(exists($domains->{$i}))){
      $islands->[$mark]=[()];
      $domains->{$i}=$mark;
      MarkMyNeighbors($R,$i,$domains,$DisplaceList,$rcut,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_);
      $mark++;
    }
  }
  for $i (keys %$domains){
    push @{$islands->[$domains->{$i}]},$i;
  }
  return $islands;
}

sub MarkMyNeighbors{ # a recursive subroutine
  my ($R,$j,$domains,$DisplaceList,$rcut,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_)=@_;
  my ($i,$m,$Rijx,$Rijy,$Rijz,$Rij);
  for $i (@$DisplaceList){
    if($i != $j){ #my enighbors
      $Rijx=$R->[$i][0] - $R->[$j][0];
      $Rijy=$R->[$i][1] - $R->[$j][1];
      $Rijz=$R->[$i][2] - $R->[$j][2];
      if($Rijx < $Ax2_){$Rijx+=$Ax;}
      if($Rijy < $Ay2_){$Rijy+=$Ay;}
      if($Rijz < $Az2_){$Rijz+=$Az;}
      if($Rijx > $Ax2){$Rijx-=$Ax;}
      if($Rijy > $Ay2){$Rijy-=$Ay;}
      if($Rijz > $Az2){$Rijz-=$Az;}
      $Rij=sqrt($Rijx*$Rijx+$Rijy*$Rijy+$Rijz*$Rijz);
      if($Rij < $rcut){
        if(!(exists($domains->{$i}))){
          $domains->{$i}=$domains->{$j};
          MarkMyNeighbors($R,$i,$domains,$DisplaceList,$rcut,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_);
        }
      }
    }# end of outer if block
  }
}

sub FindMostUndercoordianted{
  my ($i,$islands,$NCord)=@_;
  my ($mincord,$j);
  $mincord=$NCord->{$islands->[$i][0]};
  for $j (@{$islands->[$i]}){
    if($NCord->{$j} < $mincord){
      $mincord=$NCord->{$j};
    }
  }
  return $mincord;
}

sub readIslands{
  my ($islandsinfo,$rcut,$total_atoms)=@_;
  my ($good,@line,$line,$i,$j,$total,$numatoms,$islands,$mincord_island);
  open(ISMIN,"<$islandsinfo") || die "cannot open $islandsinfo\n";
  $i=$j=$good=$total=0;
  while($line=<ISMIN>){
    $line=~s/^\s+//;
    if($line eq "\n"){next;}
    @line=split(/\s+/,$line);
    if(lc($line[0]) eq "rcut"){
      if($line[1]=~/^[-+]?\d+(\.\d+)?([eE][-+]?\d+)?$/ && $line[1]==$rcut){
        $good =1;
        $i++;
        next;
      }else{
        last;
      }
    }
    if($i==1 && lc($line[0]) eq "island"){
      $j=0;
      $i++;
      next;
    }
    if($i==2){
      $numatoms=@line-1;
      $mincord_island->[$j]=$line[1];
      $islands->[$j]=[(@line[2 .. $numatoms])];
      $j++;
      $total=$total+$numatoms-1;
    }
  }
  close ISMIN;
  if($total != $total_atoms){ $good=0; }
  return ($islands,$mincord_island,$good);
}

sub writeIslands{
  my ($islandsinfo,$islands,$mincord_island,$rcut)=@_;
  my ($i,$j);
  open(ISMIN,">$islandsinfo") || die "cannot open $islandsinfo\n";
  print ISMIN "Rcut ".$rcut."\n";
  print ISMIN "Island\tMinCoord\tAtomsInIsland\n";
  for($i=0;$i<@$islands;$i++){
    print ISMIN $i."\t".$mincord_island->[$i]."\t"."@{$islands->[$i]}\n";
  }
  close ISMIN;
}
