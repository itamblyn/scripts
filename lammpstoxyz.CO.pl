#!/usr/bin/perl

###
sub nint {
  my $x = $_[0];
  my $n = int($x);
  if ( $x > 0 ) {
    if ( $x-$n > 0.5) {
      return $n+1;
    }
    else {
      return $n;
    }
  }
  else {
    if ( $n-$x > 0.5) {
      return $n-1;
    }
    else {
      return $n;
    }
  }
}
###

$atom_types = "C O";
print "WARNING!! check atom types: $atom_types\n";
@atom = split " ", $atom_types;
open (INF, "<$ARGV[0]") || die "Can't open input dump file\n";
open (OUT, ">lmp_traj.xyz") || die "Can't open lmp_traj.gen\n";
while (<INF>) {
  $time_step = <INF>;
  chomp($time_step);
  $line = <INF>;
  $natom = <INF>;
  chomp($natom);
  print (OUT "$natom\n");
  $line = <INF>;
  # assumes orthorhombic box
  @line = split " ", <INF>;
  $lx = $line[1] - $line[0];
  @line = split " ", <INF>;
  $ly = $line[1] - $line[0];
  @line = split " ", <INF>;
  $lz = $line[1] - $line[0];
  print (OUT "$lx $ly $lz # time step $time_step\n");
  @line = split " ", <INF>;
  for ($i = 0; $i < @line; $i++) {
    if ($line[$i] =~ /id/) {
      $id_ind = $i-2;
    } elsif ($line[$i] =~ /type/) {
      $type_ind = $i - 2;
    } elsif ($line[$i] =~ /x/) {
      $x_ind = $i - 2;
    } elsif ($line[$i] =~ /y/) {
      $y_ind = $i - 2;
    } elsif ($line[$i] =~ /z/) {
      $z_ind = $i - 2;
    }
  }
  for ($i = 0; $i < $natom; $i++) {
    @line = split " ", <INF>;
    $line[$x_ind] -= $lx*nint($line[$x_ind]/$lx);
    $line[$y_ind] -= $ly*nint($line[$y_ind]/$ly);
    $line[$z_ind] -= $lz*nint($line[$z_ind]/$lz);
    print (OUT "$atom[$line[$type_ind]-1] $line[$x_ind] $line[$y_ind] $line[$z_ind]\n");
  }
}
close (INF);
close (OUT);
