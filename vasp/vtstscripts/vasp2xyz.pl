eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

#-----------------------------------------------------------------------------------------
# vasp2xyz: version 1.1 (c) April 2001 by MD 
# Program to convert a POSCAR and CONTCAR files into xyz format to view it with RasMol 
#
# Help: see vasp2xyz.hlp or type `vasp2xyz -h` 
#-----------------------------------------------------------------------------------------

use FindBin qw($Bin);

$version = 1.1;

sub ReadInput {
  my @input = ("0","POSCAR","","0","0.03");
  my $ncon = 0;

  foreach (@ARGV) {
    if (/-c/) {
      $input[0] = 1;
    } 
    elsif (/-h/) {
      displayHelp ("$Bin/vasp2xyz.help");
      exit;
    } 
    elsif (/-v/) {
      print "\nvasp2xyz: $version (c) April 2001 by MD\n";
      exit;
    } 
    elsif (/-nos/) {
      $input[3] = 1;
    } 
    elsif (/^dl=/) {
      $input[4] = substr($_,3);
    } 
    elsif (/-n/) {
      $ncon = 1;
    } 
    elsif (/^i=/) {
      $input[1] = substr($_,2);
    } 
    elsif (/^o=/) {
      $input[2] = substr($_,2);
    }
  }
  if ($input[2] eq "") {
    if ($ncon == 1) {
      if ($input[1] =~ m/\./g) {
        $input[2] = substr ($input[1], pos $input[1]);
        $input[2] = $input[2].".xyz";
      } else {
        print "Error: The inputfilename is NOT compatible with -n naming option!\n";
        exit;
      }
    }
    else {
      $input[2] = $input[1].".xyz";
    }
  }
  return @input;
}

sub displayHelp ()
{
  open (HELPFILE, "<$_[0]") or die "Can't open help file: $!\n"; 
  while (<HELPFILE>) {
    print $_;
  }
  close (HELPFILE);
}

sub completeBox {
  # input list: (point,b1,b2,b3)
  my @addpos = ();
  my @newp = (0,0,0);
  my @newpdirect;
  my $dl = shift @_;
  my @basis = @_;
  my @p;
  my ($i,$j,$k);
  my @m = ();
  my @outer;
  my $inboxshell;

  $p[0] = shift @basis;
  $p[1] = shift @basis;
  $p[2] = shift @basis;

  for ($i=0;$i<2;$i++) {
    for ($j=0;$j<2;$j++) {
      for ($k=0;$k<2;$k++) {
        $newp[0] = $p[0] + $i*$basis[0] + $j*$basis[3] + $k*$basis[6];
        $newp[1] = $p[1] + $i*$basis[1] + $j*$basis[4] + $k*$basis[7];
        $newp[2] = $p[2] + $i*$basis[2] + $j*$basis[5] + $k*$basis[8];

        @newpdirect = getDirectCoord (@newp,@basis);
        @outer = inouterBox ($dl,@newpdirect);
        $inboxshell = shift @outer;

        if ($inboxshell == 1) {
          push @addpos, @newp;
	}   
      }
    }
  }    
  return @addpos;
}

sub inBox {
  # input list: (point,b1,b2,b3)
  my $in = 1;
  my @basis = @_;
  my @p = ();
  my @n = (0,0,0);
  my @m0;
  my @m1;
  my @m2;
  my $d;

  $p[0] = shift @basis;
  $p[1] = shift @basis;
  $p[2] = shift @basis;

  push @m0, @basis;
  push @m1, @basis;
  push @m2, @basis;
  
  $m0[0] = $p[0];   
  $m0[1] = $p[1];   
  $m0[2] = $p[2];   

  $m1[3] = $p[0];   
  $m1[4] = $p[1];   
  $m1[5] = $p[2];   

  $m2[6] = $p[0];   
  $m2[7] = $p[1];   
  $m2[8] = $p[2];   

  $d = det3D(@basis);

  $n[0] = det3D(@m0)/$d;   
  $n[1] = det3D(@m1)/$d;   
  $n[2] = det3D(@m2)/$d;   

  foreach (@n) {
    if ($_ > 1) {
     $in = 0; 
    }
  }
  return $in;
} 

sub getDirectCoord {
  # input list: (point,b1,b2,b3)
  my $in = 1;
  my @basis = @_;
  my @p = ();
  my @n = (0,0,0);
  my @m0;
  my @m1;
  my @m2;
  my $d;

  $p[0] = shift @basis;
  $p[1] = shift @basis;
  $p[2] = shift @basis;

  push @m0, @basis;
  push @m1, @basis;
  push @m2, @basis;
  
  $m0[0] = $p[0];   
  $m0[1] = $p[1];   
  $m0[2] = $p[2];   

  $m1[3] = $p[0];   
  $m1[4] = $p[1];   
  $m1[5] = $p[2];   

  $m2[6] = $p[0];   
  $m2[7] = $p[1];   
  $m2[8] = $p[2];   

  $d = det3D(@basis);

  if ($d == 0) {
    print "Error: Basis vectors are NOT linear independent!\n";
    die;
  } else {
    $n[0] = det3D(@m0)/$d;   
    $n[1] = det3D(@m1)/$d;   
    $n[2] = det3D(@m2)/$d;   
  }
  return @n;
} 

sub det3D {
 # 3x3 matrix is entered as the stack of three column vectors m = (a1,a2,a3)
 my $d = 0;
 my @m = @_;
 
 $d = $m[0]*($m[4]*$m[8]-$m[7]*$m[5]) + $m[3]*($m[7]*$m[2]-$m[1]*$m[8]) + $m[6]*($m[1]*$m[5]-$m[4]*$m[2]);
 return $d;
}

sub inouterBox {
  # input list: (dl, n)
  my $dl = shift @_;
  my @n = @_;
  my $inouterbox = 0;
  my ($i,$j,$k);
  my @newn = (0,0,0);

  for ($i=0; $i < 3; $i++) {
    $newn[$i] = $n[$i];
    if (($n[$i] > (1-$dl)) && ($n[$i] < (1+$dl))) {

      $j = ($i + 1) % 3;
      $k = ($i + 2) % 3;
 
      if (($n[$j] > -1*$dl) && ($n[$j] < (1+$dl)) && ($n[$k] > -1*$dl) && ($n[$k] < (1+$dl))) {
        $inouterbox = 1;
        # correct if in outer box
        $newn[$i] = $newn[$i] - 1;
      }
    }
  }
  return ($inouterbox,@newn);
}

sub insertSpecy {
  # input list: (specyname,p1,p2,p3,p4,....) 
  my $specy = shift @_;
  my @m = @_;
  my @newm = ();
  my $i = 0;

  foreach (@m) {
    if ($i % 3 == 0) {
      push @newm, $specy;
    }
    push @newm, $_;
    $i++; 
  } 
  return @newm;
}

sub scalePoints {
  # input list need to be of the form ("Si","1","1","1",......)
  my $s = shift @_;
  my @m = @_;
  my $l = scalar @m;
  
  for ($i=0;$i < $l;$i++) {
    if ($i % 4 != 0) {
      $m[$i] = $s*$m[$i];
    }
  }
  return @m;
}

# -----------------------------------------------------------------------------
# main program
# -----------------------------------------------------------------------------

@in = ReadInput ();

$infile = $in[1];
$outfile = $in[2];
$name = "";
$s = 1;
@basis = (1,0,0,0,1,0,0,0,1);
@pos = ();
@species = ();
$i=1;
$readdata = 0;
$count = 0;
$cartesiandirect = 0;

open (INFILE,"<$infile") or die "Can't open file: $!\n";

while (<INFILE>) {
  if ($i==1) {
    @species = split ' ',$_;
    $name = shift @species;
  } 
  elsif (($i==2) && ($_ > 0)) {
    $s = $_;
  }
  elsif (($i>2) && ($i<6)) {
    @b = split ' ',$_;
    $basis[($i-3)*3] = $b[0];
    $basis[($i-3)*3+1] = $b[1];
    $basis[($i-3)*3+2] = $b[2];

    if (det3D(@basis) == 0) {
      print "Error: The basis vectors in the POSCAR/CONTCAR file are NOT linear independent!\n";
      exit;
    }
  }    
  elsif ($i==6) {
    @nospecies = split ' ',$_;
    $n = 0;
    foreach (@nospecies) {
      $n = $n + $_;
    }
    $atoms = $n;
    if (scalar @nospecies != scalar @species) {
      print "Error: The number of species is NOT consistent in the POSCAR/CONTCAR file!\n";
      exit;
    }
  }
  elsif ((/^[CcKk]/) && ($readdata == 0)) {
    $cartesiandirect = 1;
    $readdata = 1;
  }
  elsif ((/^[Dd]/) && ($readdata == 0)) {
    $cartesiandirect = 2;
    $readdata = 1;
  }
  elsif ($readdata == 1) {
    @data = split ' ', $_;
    @xyz = (0,0,0);	

    # cartesian mode
    if ($cartesiandirect == 1) {
      $xyz[0] = $data[0];
      $xyz[1] = $data[1];
      $xyz[2] = $data[2]; 
    }
    # direct mode
    elsif ($cartesiandirect == 2) {
      $xyz[0] = $data[0]*$basis[0] + $data[1]*$basis[3] + $data[2]*$basis[6];
      $xyz[1] = $data[0]*$basis[1] + $data[1]*$basis[4] + $data[2]*$basis[7];
      $xyz[2] = $data[0]*$basis[2] + $data[1]*$basis[5] + $data[2]*$basis[8]; 
    }

    # fix periodic boundary box problem; $dl is acceptance range
    $dl = $in[4];
    @xyzdirect = getDirectCoord (@xyz, @basis);
    @outer = inouterBox ($dl,@xyzdirect);
    $outside = shift @outer;
    @xyzdirect = @outer;

    if ($outside == 1) {
      # print "Fix periodic boundary problem!\n";
      $xyz[0] = $xyzdirect[0]*$basis[0] + $xyzdirect[1]*$basis[3] + $xyzdirect[2]*$basis[6];
      $xyz[1] = $xyzdirect[0]*$basis[1] + $xyzdirect[1]*$basis[4] + $xyzdirect[2]*$basis[7];
      $xyz[2] = $xyzdirect[0]*$basis[2] + $xyzdirect[1]*$basis[5] + $xyzdirect[2]*$basis[8]; 
    }

    $specy = 0;
    $m = 0;
    $count++;

    if ($atoms == $count) {
      $readdata = 2;
    }

    foreach (@nospecies) {
      $m = $m + $_;
      if ($count > $m) {
        $specy++;
      }
    }

    push @pos, $species[$specy];
    push @pos, @xyz;

    # completeBox
    if ($in[0] == 1) {
      @apos = completeBox ($dl,@xyz,@basis);     

      unshift @apos, $species[$specy];      

      @newpos = insertSpecy (@apos);
      $n = $n + scalar (@newpos)/4;
      push @pos, @newpos;
     } 
  }
  $i++;
}
close (INFILE);

# Rescale all points

if ($in[3] == 0) {
  @rescale = ($s,@pos);
  @pos = ();
  push @pos, scalePoints (@rescale);
}

# Write XYZ file

open (OUTFILE, ">$outfile") or die "Can't open file: $!\n"; 
print OUTFILE "$n\n";
print OUTFILE "$name\n";
for ($i=0; $i < $n; $i++) {
  printf OUTFILE "%s\t%10.6f\t%10.6f\t%10.6f\n",$pos[$i*4],$pos[$i*4+1],$pos[$i*4+2],$pos[$i*4+3];
}
close (OUTFILE);
