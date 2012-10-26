#!/usr/bin/perl

# S.A. Bonev, LLNL, 02/13/03
#
# usage plnd.pl he i1 i2 ta 
#
# plots the number of nearest neighbour pairs the distances between which is
# larger than specified values; processes jeep output files
# ${he}${i1}${ta} to ${he}${i2}${ta};
# plots number vs. t (time), where t is the simulation time in fs;

push(@INC,"$ENV{'HOME'}/util");
require 'stat_subs.pl';
require 'gnuplot_subs.pl';

$au2fs = 0.0249;          # 1 a.u. time = 0.0249 fs
$rd1 = 1.8; $rd2 = 1.9; $rd3 = 2.0; $rd4 = 2.1;

$he = shift(@ARGV);  if ($he eq 'def') { $he = 'D128_MD.out.'; }
$i1 = shift(@ARGV);
$i2 = shift(@ARGV);
$di = 5; # $di = shift(@ARGV);
$ta = shift(@ARGV);

$datafile = 'nd.dat';

$irec = 0; @r_ave = (); @r_std = (); @t = ();     

$atom_flag=1;
for ($i = $i1; $i<=$i2; $i++) {
    open(IN, "< $he$i$ta");        # process the ouput files, one at a time 
    @X = (); @Y = (); @Z = ();

    $nsnap=0;
    while (<IN>) {
	chomp;

	if (/iprint\s+=\s+(\S+)/){ $ip = $1; }
	if (/dt\s+=\s+(\S+)\s+/) { $dt = $1; }
	if (/Total of\s+(\S+)\s+/){$N  = $1; }
	if (/run\s+(\S+)/)    { $niter = $1; } 
	if (/atom\s+(\S+)1\s+/ && $atom_flag){ $atom = $1; $atom_flag = 0;}
	if (/set_cell\s+(\S+)\s+(\S+)\s+(\S+)/) { $ca=$1; $cb=$2; $cc=$3; } 

	if (/\#\#\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/) {
	    $xx = $2; $yy = $3; $zz = $4; 
	    push @X, $xx; push @Y, $yy; push @Z, $zz; 
	    $atom_cur = $1;

	    if ($atom_cur eq "${atom}1" || $atom_cur eq "\*${atom}1" ) {   # record the snapshot time
	        $tt = $ip*$irec*$dt*$au2fs; 
	        $irec++; $nsnap++;
	        push @t, $tt; 
            }
	}
	
    }
    close(IN);

#    $nsnap = $niter/$ip;        # the output file contains nsnap snapshots

    for ($is=0; $is<$nsnap; $is=$is+$di) { 
	$nd1 = 0;  $nd2 = 0; $nd3 = 0; $nd4 = 0;
	@x=@X[$N*$is .. $N*$is+$N-1];
	@y=@Y[$N*$is .. $N*$is+$N-1];     # separate the snapshots 
	@z=@Z[$N*$is .. $N*$is+$N-1];
	for ($m=0; $m<$N; $m++) {
	    @dr = ();
	    for ($l=0; $l<$N; $l++){
		if ( $m != $l ) {         # compute all interatomic distances 
		    $dx = $x[$l]-$x[$m]; $dx = $dx - $ca*&anint($dx/$ca); 
		    $dy = $y[$l]-$y[$m]; $dy = $dy - $cb*&anint($dy/$cb);
		    $dz = $z[$l]-$z[$m]; $dz = $dz - $cc*&anint($dz/$cc);
		    push @dr, sqrt($dx*$dx + $dy*$dy + $dz*$dz);
		}
	    }
	    $dr_min = &min(@dr);      # find the n.n. distance
	    if ( $dr_min > $rd1 ) { $nd1++ ;}
	    if ( $dr_min > $rd2 ) { $nd2++ ;}
	    if ( $dr_min > $rd3 ) { $nd3++ ;}
	    if ( $dr_min > $rd4 ) { $nd4++ ;}
	}
	push @ndis1, $nd1; push @ndis2, $nd2; 
	push @ndis3, $nd3; push @ndis4, $nd4;
    }
}


open(OUT,">$datafile");
for ($i=0; $i<$irec/$di; $i++) {
    printf (OUT "%8.4f %g %g %g %g \n", $t[$di*$i], 
            $ndis1[$i], $ndis2[$i], $ndis3[$i], $ndis4[$i] );
}
close(OUT);

&plot_4($datafile,"number of nn distances > $rd1, $rd2, $rd3, $rd4 "); 



sub min {
    my $min_so_far = shift @_;
    foreach(@_) {
	if ($_ < $min_so_far) { $min_so_far = $_; }
    }
    $min_so_far;
}

sub anint {
    my $x = $_[0]; 
    my $i = int($x);
    if (abs($x-$i) >= 0.5) { 
	if ($x > 0){ $i = $i + 1; }
	if ($x < 0){ $i = $i - 1; }
    }
    $i;
}


sub plot_4 {
    my $datafile = $_[0]; my $title =  $_[1];
    open (PLOT, "|/usr/bin/gnuplot -persist");
    print PLOT <<DONE;
    set title "$title"
    plot "$datafile" using 1:2 with lines, "$datafile" using 1:3 with lines, "$datafile" using 1:4 with lines, "$datafile" using 1:5 with lines
    save "${datafile}.gnu"
DONE
    close (PLOT); 
}

