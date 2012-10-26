eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

foreach $a (@ARGV){
  if(-e "$a/OUTCAR.gz"){ system " gunzip $a/OUTCAR.gz;" ; }
  elsif("$a/OUTCAR.bz2"){ system " bunzip2 $a/OUTCAR.bz2 ;" ; }
}
