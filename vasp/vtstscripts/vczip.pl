eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

$path=`pwd`;
chop($path);

$zip = $ENV{'VTST_ZIP'} ;
if($zip eq ''){ $zip = 'gzip' ; }

foreach $a (@ARGV){
  if(-e "$a/OUTCAR") { system "$zip -9 $a/OUTCAR &" ; }
  if(-e "$a/stdout") { system "$zip -9 $a/stdout &" ; }
}
