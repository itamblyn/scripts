eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

use FindBin qw($Bin);
use lib "$Bin";
use Math::Approx;
            
sub poly {
    my($n,$x) = @_;
    return $x ** $n;
}

@args=@ARGV;
(@args>=5) || die "usage: dymanalyze.pl <flag> <displacement 1> <matrix 1> <displacement 2> <matrix 2> ...\n";

$flag=$args[0];
$order=1;
$large=5;
@args=@args[1..@args-1];

for ($i=0;$i<@args;$i=$i+2) {
    $j=0;
    
    open (IN,$args[$i+1]) || die ("Can't open file $!");
    while (<IN>) {
	$_=~s/^\s+//;
	@temp=split(/\s+/,$_);
	$kk=@temp;
	for ($k=0;$k<@temp;$k++) {
	    $matrix->[$i/2][$j][$k]=$temp[$k]*$args[$i];
	}
	$j++;
    }
    close(IN);
}
$jj=$j;

print "-------------------------------------------------------------------------------------\n";
print "Dimensions: $jj $kk    Flag: $flag\n";
if ($flag>0) {
    print "Fit $flag points, check how fitted line compares to data points (error is *100)\n";
} else {
    print "Successively fit line to more and more data points, check how force constants change\n";
}
print "-------------------------------------------------------------------------------------\n";

printf "           reference        ";
for ($i=0;$i<@args;$i=$i+2) {
    printf "%18.13f          ",$args[$i];
}
print "\n";

for ($j=0;$j<$jj;$j++) {
    for ($k=0;$k<$kk;$k++) {
	undef %x;
#	undef %y;
	$number=$flag;
	if ($number==0) {$number=1;}
	for ($i=0;$i<$number*2;$i=$i+2) {
	    $x{$args[$i]}=$matrix->[$i/2][$j][$k];
#	    print "  x: ".$args[$i]." ".$x{$args[$i]}." -- ";
	}
	$x{0}=0;
	$fit = new Math::Approx (\&poly, $order, %x);

	$fit1=${$fit->{'A'}}[1];
    
    if (abs($fit1)>$large) {
#	$fit->print;
	printf  "%3i %3i %18.13f  ",($j+1),($k+1),$fit1."  ";
    }

    undef %y;
    for ($i=0;$i<@args;$i=$i+2) {
	if ($flag==0) {
#	    undef %y;
	    $y{0}=0;
	    $y{$args[$i]}=$matrix->[$i/2][$j][$k];
	    $fit=new Math::Approx (\&poly, $order, %y);
	    $fit2=${$fit->{'A'}}[1];
#	if (abs($fit1)>$large) {
#	    $fit->print;
#	}
	$error=abs($fit1-$fit2)/abs($fit1)*100;
#	print "a";
    } else {
	$fit2=$fit->approx($args[$i]);
#	$error=abs($matrix->[$i/2][$j][$k]-$fit2)/abs($fit2)*100;
	$error=abs($matrix->[$i/2][$j][$k]-$fit2)*100;
#	print "b";
    }
    if (abs($fit1)>$large) {
	printf "%18.13f %6.3f\%  ",$fit2,$error;
	$errornum->[$i/2]++;
	$errortrack->[$i/2][$errornum->[$i/2]]=$error;
	$errorsum->[$i/2]+=$error;
#	$errornum->[$i/2]++;
    }
}
#if (abs($fit1)>$large) {
#    $fit->print;
#}
	
#	printf  "%18.13f", ${$fit->{'A'}}[1]."  ";
if (abs($fit1)>$large) {
    print "\n";
}
    }
#    print "\n";
}

for ($i=0;$i<@args;$i=$i+2) {
    for ($j=1;$j<=$errornum->[$i/2];$j++) {
    $variance->[$i/2]+=($errorsum->[$i/2]/$errornum->[$i/2]-$errortrack->[$i/2][$j])**2;
}
$variance->[$i/2]=sqrt($variance->[$i/2]/($errornum->[$i/2]-1));
}

printf "Flag: %2i: Disp, ave error, sd: ",$flag;
for ($i=0;$i<@args;$i=$i+2) {
    printf "%5.3f %6.3f\% ",$args[$i],($errorsum->[$i/2]/$errornum->[$i/2]);
printf "(%6.3f) -- ",$variance->[$i/2];
}
print "\n";
