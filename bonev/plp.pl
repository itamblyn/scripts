#!/usr/bin/perl

# S.A. Bonev, LLNL, 01/21/03
#
# usage plp.pl he i1 i2 ta 
#
# plots a running average of  pressure from jeep output files
# ${he}${i1}${ta} to ${he}${i2}${ta};
# plots it in P vs. t (time), where P is in GPa and t is in fs;

push(@INC,"$ENV{'HOME'}/util");
require 'stat_subs.pl';
require 'gnuplot_subs.pl';

$au2fs = 0.0249;  # 1 a.u. time = 0.0249 fs
$t_av  = 40;      # av_t in fs; P(t) is averaged from t-0.5*t_av to t+0.5*t_av


$he = shift(@ARGV);  if ($he eq 'def') { $he = 'D128_MD.out.'; } 
$i1 = shift(@ARGV);
$i2 = shift(@ARGV);
$ta = shift(@ARGV);

$datafile = "P_ra_${t_av}.dat";

$irec = 0; @P = (); @t = ();

for ($i = $i1; $i<=$i2; $i++) {
    open(IN, "< $he$i$ta");
    while (<IN>) {
	chomp;
	if (/iprint\s+=\s+(\S+)/){ $ip = $1; }
	if (/dt\s+=\s+(\S+)\s+/) { $dt = $1; }
	if (/Pressure\s+\(GPa\)\s+(\S+)/) {
	    $pp = $1; $tt = $ip*$irec*$dt*$au2fs;
	    $irec++;
	    push @P, $pp; push @t, $tt;
	}
    }
    close(IN);
}

$i_av = int($t_av/(2*$dt*$ip*$au2fs));

@P_ra = ();
open(OUT,">$datafile");
for ($i = $i_av; $i <= $irec-$i_av; $i++) {
    @P_seg = @P[$i-$i_av .. $i+$i_av];
    $pp  = &average(@P_seg);  push @P_ra, $pp;
}

$i_ra = $irec - 2*$i_av;
open(OUT,">$datafile");
for ($i=0; $i<$i_ra; $i++) {
    printf (OUT  "%4.0f %6.3f\n", $t[$i]+$t[$i_av], $P_ra[$i]); 
}
close(OUT);

$P_ave = sprintf "%.2f", &average(@P_ra);
$P_std = sprintf "%.2f", &stdev(@P_ra);


&plot_simple_w_title($datafile,"P = ($P_ave +/- $P_std) GPa");


