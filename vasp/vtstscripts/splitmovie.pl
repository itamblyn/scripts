eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

# 07-05-2004

# Because Jmol is so stupid it is necessary to brake the movie.xyz file
# down to smallar pieces sometimes.

# Check if there is a movie file there, and if so if it is old
  if(-e "movie.xyz"){
    if((-M "XDATCAR") < (-M "movie.xyz")){
#     It is time to generate a new movie file 
      system "xdat2xyz.pl" ;    
    } 
  }else{
#   There is NO movie file
    system "xdat2xyz.pl" ; 
  }

# Get the number of atoms and the length of the xyz file
  $natoms = `head -1 movie.xyz` ;
  chomp($natoms) ;
  $nl = `wc -l movie.xyz` ;
  $nl=~s/^\s+//g ; @nl=split /\s+/,$nl ; 
  $nl = @nl[0] ;
  $nsteps = $nl/($natoms+2) ;

  print "\n" ;
  print "  The number of atoms is  :  ",$natoms,"\n" ;
  print "  The number of steps is :   ",$nsteps,"\n" ;
  print "\n" ;
  print "  Into how many pieces should the movie file been split upto? ... " ;
  $np = <STDIN> ;  
  $base = int($nsteps/$np) ;
  print "  The last file will have ",$base+($nsteps-$base*$np)," steps","\n" ;
  print "  All other files will have ",$base," steps","\n" ;
  for($i=0; $i<$np-1; $i++){
    $spl .= $base*($i+1)*($natoms+2)+1 ;
    $spl .= " " ;
  }
  system "csplit movie.xyz $spl >& /dev/null" ;
