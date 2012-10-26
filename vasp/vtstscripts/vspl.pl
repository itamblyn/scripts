eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#;-*- Perl -*-

@args=@ARGV;

# default values for the aguments

$dJ=20;
$NumI=scalar(@args);
$NumJ=($NumI-1)*$dJ+1;

print "$NumI\n";
print "Making movie with $NumJ images.\n";

# interpolate the NumJ images

system("mkdir mov"); 

$j=0;
for($i=1; $i<$NumI; $i++){
        print "Image ",$i,"\n";
	$file1=$args[$i-1];
	$file2=$args[$i];
        print "File1 ",$file1,"\n";
        print "File2 ",$file2,"\n";
	$command="vasp2con.pl $file1;";
	$command.="vasp2con.pl POSCAR.con ciPOSCAR;";
	$command.="mv ciPOSCAR p1;";
	$command.="vasp2con.pl $file2;";
	$command.="vasp2con.pl POSCAR.con ciPOSCAR;";
	$command.="mv ciPOSCAR p2;";
	system("$command");
	for($dj=0; $dj<$dJ; $dj++){
		$fract=$dj/$dJ;
                print "Frac ",$fract,"\n";
		$command="interpolate.pl p1 p2 $fract;";
#		$command.="quad.pl POSCAR.out POSCAR;";
		$command.="mv POSCAR.out POSCAR;";
		$command.="vasp2con.pl POSCAR;";
		$command.="mv POSCAR.con mov/";
		if($j<10){$filename="img000"."$j".".con";}
		elsif($j<100){$filename="img00"."$j".".con";}
		elsif($j<1000){$filename="img0"."$j".".con";}
		else{$filename="img"."$j".".con";}
		$command.="$filename;";
		system("$command");
		$j++;
	}
}

$file1=$args[$NumI-1];
$command="vasp2con.pl $file1;";
$command.="vasp2con.pl POSCAR.con ciPOSCAR;";
$command.="vasp2con.pl ciPOSCAR;";
$command.="mv POSCAR.con mov/";
if($j<10){$filename="img000"."$j".".con";}
elsif($j<100){$filename="img00"."$j".".con";}
elsif($j<1000){$filename="img0"."$j".".con";}
else{$filename="img"."$j".".con";}
$command.="$filename;";
$command.="rm ciPOSCAR;";
$command.="rm POSCAR; rm p1; rm p2;";
system("$command");
