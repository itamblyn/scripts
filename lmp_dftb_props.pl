#!/usr/bin/perl
use constant atm2gpa => 0.000101325; 

open (INF, "<$ARGV[0]") || die "Can't open input: $ARGV[0]\n";
open (OUT, ">lmp.log")  || die "Can't open lmp.log\n";
open (SCC,  ">scc.dat") || die "Can't open scc.dat\n";
open (ETOT, ">etot_lmp.out") || die "Can't open etot_lmp.out\n";

if ($ARGV[1] =~ /^\d+$/) {
  $cutoff = $ARGV[1];
  print "Averaging with cutoff: $cutoff\n"; 
}
if ($ARGV[2]) {
  $stop_time = $ARGV[2];
} 
$scc_sample = "8   -0.32231779E+03   -0.20182474E-02    0.13764296E+00\n";

$avg_ptot = 0;
$avg_pxx = 0;
$avg_pyy = 0;
$avg_pzz = 0;
$avg_lx = 0;
$avg_ly = 0;
$avg_lz = 0;
$avg_temp = 0;
$avg_ener = 0;
$avg_pe = 0;
$avg_ke = 0;
$avg_dhug = 0;
$avg_dray = 0;
$step_counter = 0;
$cutoff_count = 0;
$lx = 0;
$first = 1;
while (<INF>) {
  $etot = 0;
  @line = split " ";
  $len = @line;
  if (/atoms/) {
    @line = split " ";
    if ($line[0] =~ /^\d+$/) {
      $natoms = $line[0];
    }
  }
  if ((/Step/) && ($first)) {
    $dat_length = $len;
    $first = 0;
    for ($i = 0; $i < $len; $i++) {
      if ($line[$i] =~/Step/) {
        $step_ind = $i;
      } elsif ($line[$i] =~ /Temp/) {
        open (TEMP, ">temperature.dat") || die "Can't open temperature.dat\n";
        $temp_ind = $i;
      } elsif ($line[$i] =~ /PotEng/) {
        open (PE,   ">pe.dat")  || die "Can't open pe.dat\n";
        $pe_ind = $i;
      } elsif ($line[$i] =~ /KinEng/) {
        open (KE,   ">ke.dat")  || die "Can't open ke.dat\n";
        $ke_ind = $i;
      } elsif ($line[$i] =~ /TotEng/) {
        open (ENER, ">ener.dat") || die "Can't open ener.dat\n";
        $ener_ind = $i;
      } elsif ($line[$i] =~ /Press/) {
        open (PTOT, ">ptot.dat") || die "Can't open ptot.dat\n";
        $press_ind = $i;
      } elsif ($line[$i] =~ /Pxx/) {
        open (PXX,  ">pxx.dat") || die "Can't open pxx.dat\n";
        $pxx_ind = $i;
      } elsif ($line[$i] =~ /Pyy/) {
        open (PYY,  ">pyy.dat") || die "Can't open pyy.dat\n";
        $pyy_ind = $i;
      } elsif ($line[$i] =~ /Pzz/) {
        open (PZZ,  ">pzz.dat") || die "Can't open pzz.dat\n";
        $pzz_ind = $i;
      } elsif ($line[$i] =~ /Lx/) {
        open (LX,   ">lx.dat")  || die "Can't open lx.dat\n";
        $lx_ind = $i;
      } elsif ($line[$i] =~ /Ly/) {
        $ly_ind = $i;
      } elsif ($line[$i] =~ /Lz/) {
        $lz_ind = $i;
      } elsif ($line[$i] =~ /lgr_vel/) {
        open (LVEL, ">lgr_vel.dat") || die "Can't open lgr_vel.dat\n";
        $lgr_vel_ind = $i;
      } elsif ($line[$i] =~ /dhug/) {
        open (DHUG, ">dhug.dat") || die "Can't open dhug.dat\n";
        $dhug_ind = $i;
      } elsif ($line[$i] =~ /dray/) {
        open (DRAY, ">dray.dat") || die "Can't open dray.dat\n";
        $dray_ind = $i;
      }
    }
  }
  if (($len == $dat_length) && ($line[0] =~ /^\d+$/) && ($line[1] =~ /\d/)) {
    $step_count = $line[0];
    $step_counter++;
    if ($step_count == $cutoff) {
      $cutoff_counter = $step_counter;
    }
    print (OUT "@line\n");
    if ($stop_time) {
      if ($step_count > $stop_time) {
        print "STOPPING ANALYSIS: $step_count\n";
        last;
      }
    }
    if ($temp_ind) {
      print (TEMP "$line[$step_ind] $line[$temp_ind]\n");
      if ($step_count > $cutoff) {
        $avg_temp += $line[$temp_ind];
      }
    }
    if ($pe_ind) {
      print (PE "$line[$step_ind] $line[$pe_ind]\n");
      if ($step_count > $cutoff) {
        $avg_pe += $line[$pe_ind];
      }
    }
    if ($ke_ind) {
      print (KE "$line[$step_ind] $line[$ke_ind]\n");
      if ($step_count > $cutoff) {
        $avg_ke += $line[$ke_ind];
      }
    }
    if ( ($ke_ind) && ($pe_ind) ) {
      $etot = $line[$pe_ind] + $line[$ke_ind];
      print (ETOT "$line[$step_ind] $etot\n");
    }
    if ($ener_ind) {
      print (ENER "$line[$step_ind] $line[$ener_ind]\n");
      if ($step_count > $cutoff) {
        $avg_ener += $line[$ener_ind];
      }
    }
    if ($press_ind) {
      $line[$press_ind] *= atm2gpa;
      if ($step_count > $cutoff) {
        $avg_ptot += $line[$press_ind]; 
      }
      print (PTOT "$line[$step_ind] $line[$press_ind]\n");
    }
    if ($pxx_ind) {
      $line[$pxx_ind] *= atm2gpa;
      if ($step_count > $cutoff) {
        $avg_pxx += $line[$pxx_ind];
      }
      print (PXX "$line[$step_ind] $line[$pxx_ind]\n");
    }
    if ($pyy_ind) {
      $line[$pyy_ind] *= atm2gpa;
      if ($step_count > $cutoff) {
        $avg_pyy += $line[$pyy_ind];
      }
      print (PYY "$line[$step_ind] $line[$pyy_ind]\n");
    }
    if ($pzz_ind) {
      $line[$pzz_ind] *= atm2gpa;
      if ($step_count > $cutoff) {
        $avg_pzz += $line[$pzz_ind];
      }
      print (PZZ "$line[$step_ind] $line[$pzz_ind]\n");
    }
    if ($lx_ind) {
      if ($step_count > $cutoff) {
        $avg_lx += $line[$lx_ind];
        $avg_ly += $line[$ly_ind];
        $avg_lz += $line[$lz_ind];
      }
      if ($step_counter == 1) {
        $lx0 = $line[$lx_ind];
        $ly0 = $line[$ly_ind];
        $lz0 = $line[$lz_ind];
      }
      print (LX "$line[$step_ind] $line[$lx_ind] $line[$ly_ind] $line[$lz_ind]\n");
    }
    if ($lgr_vel_ind) {
      print (LVEL "$line[$step_ind] $line[$lgr_vel_ind]\n");
    }
    if ($dhug_ind) {
      if ($step_count > $cutoff) {
        $avg_dhug += $line[$dhug_ind];
      }
      print (DHUG "$line[$step_ind] $line[$dhug_ind]\n");
    }
    if ($dray_ind) {
      if ($step_count > $cutoff) {
        $avg_dray += $line[$dray_ind];
      }
      print (DRAY "$line[$step_ind] $line[$dray_ind]\n");
    }
  }
  if (/iSCC Total electronic   Diff electronic      SCC error/) {
    while (<INF>) {
      if (/^\s+$/) {
        last;
      }
      $scc_data = $_;
    }
    @line = split " ", $scc_data;
    print (SCC "@line\n");
    $line = <INF>;
  }
}
close (INF);
close (OUT);
close (TEMP);
close (ENER);
close (PTOT);
close (PXX);
close (PYY);
close (PZZ);
close (LX);
$avg_ptot /= ($step_counter-$cutoff_counter);
$avg_pxx /= ($step_counter-$cutoff_counter);
$avg_pyy /= ($step_counter-$cutoff_counter);
$avg_pzz /= ($step_counter-$cutoff_counter);
$avg_temp /= ($step_counter-$cutoff_counter);
$avg_pe /= ($step_counter-$cutoff_counter);
$avg_ener /= ($step_counter-$cutoff_counter);
$avg_ke /= ($step_counter-$cutoff_counter);
$avg_lx /= ($step_counter-$cutoff_counter);
$avg_ly /= ($step_counter-$cutoff_counter);
$avg_lz /= ($step_counter-$cutoff_counter);
$avg_dhug /= ($step_counter-$cutoff_counter);
$avg_dray /= ($step_counter-$cutoff_counter);
if ($avg_ptot != 0) {
  print "P_avg (GPa): $avg_ptot\n"; 
}
print "Pxx_avg (GPa): $avg_pxx\n"; 
print "Pyy_avg (GPa): $avg_pyy\n"; 
print "Pzz_avg (GPa): $avg_pzz\n"; 
if ($avg_lx > 0) {
  print "Lx_avg (AA): $avg_lx\n"; 
  $ratio_x = $avg_lx/$lx0;
}
if ($avg_ly > 0) {
  print "Ly_avg (AA): $avg_ly\n"; 
  $ratio_y = $avg_ly/$ly0;
}
if ($avg_lz > 0) {
  print "Lz_avg (AA): $avg_lz\n"; 
  $ratio_z = $avg_lz/$lz0;
}
$ratio = $ratio_x*$ratio_y*$ratio_z;
print "V/V0 = $ratio\n";
if ($avg_dhug != 0) {
  print "DHug (K): $avg_dhug\n"; 
}
if ($avg_dray != 0) {
  print "DRay (atm): $avg_dray\n"; 
}
print "T_avg  (K):   $avg_temp\n"; 
print "PE_avg (kcal/mol):  $avg_pe\n"; 
print "KE_avg (kcal/mol):  $avg_ke\n"; 
print "Ener_avg (kcal/mol):  $avg_ener\n"; 
$avg_ener /= 23.060;
$avg_ener /= $natoms;
print "Ener_avg (eV/atom):  $avg_ener\n"; 
print "Total number of simulation time steps: $step_count\n";
$avg_step = $step_counter-$cutoff_counter;
print "Total number of time steps for averaging: $avg_step\n";
