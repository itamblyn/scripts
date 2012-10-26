eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

use FindBin qw($Bin);

system "grep ut insout.dat | grep -v itr | cut -c 5-100 > o.u.t.t.e.m.p";
system "gnuplot $Bin/insplot.gnu";
system "rm o.u.t.t.e.m.p";
