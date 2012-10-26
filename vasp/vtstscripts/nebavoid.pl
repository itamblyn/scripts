eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-
use FindBin qw($Bin);
use lib "$Bin";
use Vasp;

################################################
# Generating neb images by interpolation can create images in which two (or more) atoms 
# are too close to each other. This script applies a repelling force between each pair
# and pushes them apart until the minimum distance is reached. 
# However, this might disturb the equal space between neighboring images. We may fix it later.
# Written by Lijun Xu and Graeme Henkelman in the Univeristy of Texas at Austin on March 6, 2007
# Last modified by Lijun on March 8, 2007
################################################

if(@ARGV>0) {$rcut=$ARGV[0]; print "The minimum distance between two atoms is $rcut Angstrom.\n";}
else{die "Usage: nebavoid.pl distance (the minimum distance between two atoms; recommended value is 0.1 Angstrom.\n";}
chomp($folder=`ls -d [0-9][0-9]`);
$folder=~s/^\s+//;
@folders=split(/\s+/,$folder);
$num=@folders-1;
for($i=1;$i<$num;$i++){
  $targetfolder=$folders[$i];
  $dummy=substr($targetfolder, -1 , 1);
  if($dummy eq "\/") {chop($targetfolder);}
  if(!(-e $targetfolder)){die "folder $targetfolder does not exist.\n";;}
  $POSCAR=$targetfolder."/POSCAR";
  if(!(-e $POSCAR)) {die "file $POSCAR does not exist.\n";}
  print "working on $targetfolder...\n";
  ($coordinates,$basis,$lattice,$num_atom,$total_atoms,$selectiveflag,$selective,$description)=read_poscar($POSCAR);
  spring_relaxation($coordinates,$basis,$lattice,$total_atoms,$selective,$rcut);
  system "cd $targetfolder; mv POSCAR POSCAR_orig";
  write_poscar($coordinates,$basis,$lattice,$num_atom,$total_atoms,$selectiveflag,$selective,$description,$POSCAR);
}
sub spring_relaxation{
  my ($R,$basis,$lattice,$totalatoms,$selective,$cutoff)=@_;
  my $epsilon=0.1;
  my $stepmax=10000;
  my $drmax=0.1;
  my @line=();
  my ($i,$j,$Rijx,$Rijy,$Rijz,$Ax,$Ay,$Az,$Ax2,$Ay2,$Az2,$Ax2_,$Ay2_,$Az2_,$Rij);
  my ($force,$vector1,$vector2,$vector3,$fij,$phiij,$step,$goodsound,$converged,$frozen,$dr);
  # direct to cartesian
  $R=dirkar($R,$basis,$lattice,$totalatoms);  
  # figure out box dimension
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
    print "Fully relaxed after $step (1: trivial) steps\n";
  }else{
    print "Not fully relaxed after $stepmax steps, but we can't wait, so move on.\n";
  }
  # Cartesian to direct
  $R=kardir($R,$basis,$lattice,$totalatoms);
}

sub linear_repulsion_pot{ # E=kx-b, k< 0
  my ($r,$cutoff)=@_;
  my $slope=-0.1;
  my ($energy, $fij);
  my $intercept=$slope*$cutoff;
  $energy=$slope*$r-$intercept;
  $fij=$slope;
  return ($fij,$energy);
} 

