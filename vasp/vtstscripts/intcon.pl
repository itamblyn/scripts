eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

# This program interpolates between two con files by the given fraction.

@args=@ARGV;
@args==3 || die "usage: intcon.pl <confile1> <confile2> <fractional distance between>\n";

$confile1=$args[0];
$confile2=$args[1];
$fraction=$args[2];

$command="vasp2con.pl $confile1 ciPOSCAR;";
$command.="mv ciPOSCAR p1;";
$command.="vasp2con.pl $confile2 ciPOSCAR;";
$command.="mv ciPOSCAR p2;";
$command.="interpolate.pl p1 p2 $fraction;";
$command.="rm p1;";
$command.="rm p2;";
$command.="vasp2con.pl POSCAR.out POSCAR.con;";
$command.="rm POSCAR.out;";
system("$command");

