eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

  die "input a POSCAR or a CONTCAR \n" if @ARGV > 1 ;
  $pos = $ARGV[0] ;  
  open IN , $pos ;
  open OUT , ">$pos.vasp" ;

# Read in the input file
  while (<IN>){$file .= $_ ;}
  close IN ;  
  @file = split /\n/ , $file ;

# Get the number of names of atomi types
  $line = $file[0] ;
  chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
  $types = join " " , @line ;
  $nel = @line ;
  $line = $file[5] ;
  chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
  $numbers = join "  " , @line ;
# Calculate the total number of atoms
  while($line[$k] != undef){$natoms+=$line[$k++] ;}
  
# Write out the .vasp file
  print OUT "SOME CRAP WITH DONT CARE ABOUT   ","NCLASS=",$nel,"  ATOM=",$types,"\n" ;
  print OUT "   ",$numbers,"\n" ;
  print OUT "Direct","\n" ;
  print OUT "   ","\n" ;
  for ($i=1; $i<5; $i++){print OUT $file[$i],"\n" ;}
  print OUT "   ","\n" ;
  $sh = 7 ;
  for ($i=1;$i<=$natoms; $i++){
    $line = $file[$sh+$i] ; chomp($line) ; $line=~s/^\s+//g ; @line=split /\s+/,$line ;
    printf OUT "%13.8f %11.8f %11.8f %5s\n",@line[0..2],"#$i" ;
  }

  close OUT ;

