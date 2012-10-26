eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

foreach $a (@ARGV) {
	$command="~/bin/con2xyz.pl $a;";
	system("$command");
}
$command="cat i*.xyz > mov.xyz;";
#$command.="mv mov.xyz ..; cd ..; rm -r mov;";
system("$command");
