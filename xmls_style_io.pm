package xmls_style_io;
#use strict;
#use warnings;
use IO::File;

sub new{
  my $class = shift;
  $self = {};

  $self->{fh} = shift;

  bless $self, $class;
}

sub print_head{
  my $self = shift;
  my $fh = $self->{fh};
  $fh->print("<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<qbox:species xmlns:qbox=\"http://www.llnl.gov/casc/fpmd/qbox/ns/qbox-1.0\" 
  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" 
  xsi:schemaLocation=\"http://www.llnl.gov/casc/fpmd/qbox/ns/qbox-1.0 
  species.xsd\">\n");
}

sub xmlprint{
  my $self = shift;
  my $type = shift;
  my @content = @_;
  my $fh = $self->{fh};
  $fh->print("<${type}>\n");
  foreach $content (@content){
      $fh->print(" $content\n");
  }
  $fh->print("</${type}>\n");
}

sub xmlprintf{
  my $self = shift;
  my $type = shift;
  my @content = @_;
  my $fh = $self->{fh};
  $fh->print("<${type}>\n");
  foreach $content (@content){
      $fh->printf(" %.8f\n",$content);
  }
  $fh->print("</${type}>\n");
}

sub open_brace{
  my $self = shift;
  my $type = shift;
  my $fh = $self->{fh};
  $fh->print("<${type}>\n");
}

sub close_brace{
  my $self = shift;
  my $type = shift;
  my $fh = $self->{fh};
  $fh->print("</${type}>\n");
}

1;
