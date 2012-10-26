eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-
###################################################################
# this file was to pull poscar info from OUTCAR and place in files 
###################################################################
$var[0]=0;
$var[1]=0;
$var[2]=0;
@args=@ARGV;
$inputfilename=$args[0];
if (@args==4) {
  $var[0]=$args[1];
  $var[1]=$args[2];
  $var[2]=$args[3];
};
print "shift: @var \n";
###################################################################
  open(VASP,"POSCAR") || die "cannot open POSCAR for input \n";
  @lines=<VASP>;
  close(VASP);
  foreach $line (@lines) {chop ($line)};
  @kinds=split(' ',$lines[0]);
################################################################################
#   get lattice constant and vectors
################################################################################
  $l_con=$lines[1];
  for ($i=0;$i<3;$i++) {
    $line=$lines[$i+2];
    $line=~s/^\s+//;
    @line=split(/\s+/,$line);
    for ($j=0;$j<3;$j++) {
      $l_vec[$i][$j]=$line[$j]*$l_con;
      $d_vec[$i]+=$l_vec[$i][$j]**2;
    };
    $d_vec[$i]=sqrt($d_vec[$i]);
  };
  if ($d_vec[0]==$d_vec[1] && $d_vec[1]==$d_vec[2]) {$cubic=1};
################################################################################
#   get number of atoms
################################################################################
  $n_kinds=$lines[5];
  $n_kinds=~s/^\s+//;
  @n_kinds=split(/\s+/,$n_kinds);
  if (@n_kinds!=@kinds) {
    die "# of atom labels is not equal to number of atoms in $infile\n";
  } else {
    $n_types=@kinds;
    $count=0;
    for ($i=0;$i<@n_kinds;$i++) {
      for ($j=0;$j<$n_kinds[$i];$j++) {
        $types[$count]=$kinds[$i];
        $n_types[$count]='t'.$i;
        $count+=1;
      };
    };
  };
  $n_atoms=0;        
  for ($i=0;$i<@n_kinds;$i++) {
    $n_atoms+=$n_kinds[$i];
  };
  $inde=8;
  for ($i=$inde;$i<$inde+$n_atoms;$i++) {
    @line=split(/\b\s+/,$lines[$i]);
    $selective[$i-$inde]=$line[3]." ".$line[4]." ".$line[5];
  };
################################################################################
# print header to movie file:
################################################################################
  $outfile=">movie.vasp";
  open(VASPo,">$outfile") || die "error writing to outfile \n";
  print VASPo "Czech_yore_date_a_Fuerst   NCLASS=$n_types   ATOM=@kinds\n";
  print VASPo " @n_kinds\n";
  print VASPo "Direct\n";
################################################################################
  open(DROGAR,"OUTCAR") || die "cannot open OUTCAR for input \n";
  @lines=<DROGAR>;
  close(DROGAR);
  foreach $line (@lines) {chop ($line)};
###################################################################
  $n_step=0;
  $n_oflines=@lines;
  for($xnxn=0;$xnxn<$n_oflines;$xnxn++){
    $_=$lines[$xnxn];
    chop($_);
    @line = split;
    if ($line[0] eq "POSITION"){
      if ($line[1] eq "TOTAL-FORCE"){
        for ($j_babi=0;$j_babi<$n_atoms;$j_babi++) {
          $index=$xnxn + $j_babi +2 ;
          $_=$lines[$index];
          chop($_);
          @xline = split;
          for ($j=0;$j<3;$j++) {
            $cc[$j_babi][$j]=$xline[$j];
          };
        }
   #   print "$cc[1][1] $cc[1][2] $cc[1][3] $selective[1] \n";
#################################################################################
##   set up lattice and direct coordinates
#################################################################################
        $l_con = 1.0;
        for ($j=0;$j<$n_atoms;$j++) {
          for ($k=0;$k<3;$k++) {
            $cc_kind[$j][$k]=(($cc[$j][$k])/$l_vec[$k][$k])+$var[$k];
            if ($cc_kind[$j][$k]>1) {$cc_kind[$j][$k]--};
          };
        };
#################################################################################
##   print step image to movie file:
################################################################################

        printf VASPo "\n %.17f\n", $l_con;
        for ($i=0;$i<3;$i++) {
          printf VASPo " %22.16f%22.16f%22.16f\n",@{$l_vec[$i]}};
        printf VASPo "\n";

        for ($i=0;$i<$n_atoms;$i++) {
          printf VASPo "%20.16f%20.16f%20.16f\n", @{$cc_kind[$i]};
        };
###################################################################
      }
    }
  }
###################################################################
close(VASPo);

