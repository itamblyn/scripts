eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

# 23-01-2003

# Change the first line in POSCAR so the element symbols
# are there in the same order as in POTCAR.

  open POS, "POSCAR" 
    or die "No POSCAR to open\n" ;
  $inline = join " ", @ARGV ;
  while (<POS>) {$file.=$_;}
  close POS ;
  @file = split /\n/ , $file ;
  @file[0] = $inline,"\n" ;
  $file = join "\n" , @file ;
  open POS , ">POSCAR" ;
  print POS $file,"\n" ;

# Check to see if there is a CONTCAR there. If one is found
# change it ads well.

  if (-e "CONTCAR"){
    open POS, "CONTCAR" ;
    while (<POS>) {$file.=$_;}
    close POS ;
    @file = split /\n/ , $file ;
    @file[0] = $inline,"\n" ;
    $file = join "\n" , @file ;
    open POS , ">CONTCAR" ;
    print POS $file,"\n" ;
  } 
