#!/usr/bin/perl
package main;
$inc = "$ENV{'HOME'}/bin";
#$inc = "$ENV{'HOME'}/qbox";
push(@INC,$inc);
require 'jeep_ps.pm';
require 'xmls_style_io.pm';

#use strict;
#use warnings;
use IO::File;

die qq{Usage: ps_jeep_to_qbox.pl filename\n} unless @ARGV;
$file = shift;

$ps = new jeep_ps($file);

$out = "$file.xml";
$fh = new IO::File;
$fh->open(">$out") or die qq{Could not open $out\n};
$qbox = new xmls_style_io($fh);

$qbox->print_head();

$description = $file;
$description =~ s/_/ /g; 
$qbox->xmlprint("description",$description);

@tmp = split(/\s+/,$description);
$elem = shift(@tmp);
$qbox->xmlprint("symbol",$elem);

print "type atomic number: ";
chomp($atomic_number = <STDIN>);
$qbox->xmlprint("atomic_number",$atomic_number);

$qbox->xmlprint("mass",$ps->get_mass());

$qbox->open_brace("norm_conserving_pseudopotential");

$qbox->xmlprint("valence_charge",int($ps->get_num_valence()));

$qbox->xmlprint("lmax",$ps->get_max_l());

$qbox->xmlprint("llocal",$ps->get_l_loc());

$qbox->xmlprint("nquad",$ps->get_nquad());

if ($ps->get_nquad() == 0){
  $rquad = 0.0;
}else{
  print "type rquad: ";
  chop($rquad = <STDIN>);
}

$qbox->xmlprint("rquad",$rquad);

$qbox->xmlprint("mesh_spacing",$ps->get_mesh_spacing());

for($l = 0; $l <= $ps->get_max_l(); $l++){
  $mesh = $ps->get_num_grid($l);
  $brace = "projector l=\"$l\" size=\"$mesh\"";
  $qbox->open_brace($brace);

  @pot = $ps->get_ps($l);
#  $pot = join("\n",@pot);
#  $pot = "";
#  foreach(@pot){
#      $pot .= " " . $_ . "\n";
#  }
  $qbox->xmlprint("radial_potential",@pot);

  @wf = $ps->get_wavefun($l);
#  $wf = "";
#  foreach(@wf){
#      $wf .= " " . $_ . "\n";
#  }
#  $wf = join("\n",@wf);
  $qbox->xmlprint("radial_function",@wf);

  $qbox->close_brace("projector");
}
$qbox->close_brace("norm_conserving_pseudopotential");
$qbox->close_brace("qbox:species");

$fh->close();


