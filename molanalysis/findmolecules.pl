#! /usr/bin/perl 
#
# A program to find all the concentration of molecules in a dftb
# simulation.  First run, molanal.pl to analyze all of the files.
# Then, findmolecules.pl will find the each species present.
#
# Usage: findmolecules.pl <molanal.out>
#
# <molanal.out> is an output file from molanal.new.c
#
##################################
#  USER-DEFINED VARIABLES
#
# The block_size is the number of steps to average concentrations over.
$block_size = 100000 ;
# Molecules over this size will be counted as a "polymer".
$max_molecule_size = 48 ;
# Stop reading after this many frames. Or when end of file is found.
$stop_frame = 100000 ;
# Do not print out molecules with concentration less than this (in units of mol/cc).
$conc_floor = 1.0e-04 ;
# Minimum lifetime for a molecule (in ps).
$min_mol_life = 0.0;
# Number of iterations to use in finding reactions.
$react_iter_max = 10 ;
# Set to 1 in order to analyze rates, 0 otherwise.
$analyze_rates = 1 ;
# Reaction histories will be printed only for reactions with more than 
# this number of reaction events.
$reaction_count_floor = 10 ;
# If set to 1, calculate reverse reaction pairs.  This is somewhat slow right now.
$find_reverse_rxn = 0 ;
# Maximum charge used in calculating the histogram.
$qmax = 4.0 ;
# Charge histogram spacing.
$dq = 0.05 ;
# If 1, print out a histogram of species lifetimes.
$print_lifetime_histogram = 0 ;
#  END USER-DEFINED VARIABLES
##################################

## Turn off undeclared global variables.
use strict ;

## List of global variables.
use vars 
( 
  '$analyze_rates',       # If 1, analyze reaction rates.
  '$block_size',          # Number of steps used in block averaging over time.
  '$block_count',         # The current time block used for averaging.
  '$block_count1',        # The current time block used for averaging reactions. Necessary for consistent tracking
                          # of concentrations of reactants (which come from the step before the current one and
                          # reactions (which are detected on the current step).
  '$changed_left',        # If 1, the left hand side of the reaction changed.
  '$changed_right',       # If 1, the right hand side of the reaction changed.
  '@charge_histo',            # Charge histogram file.
  '$conc_floor',          # Only print species above this concentration (in mol/cc).
  '$find_reverse_rxn',    # If 1, find reverse reaction rates.
  '$found_charge',        # If 1, a charge file was found. 
  '$found_timestep',      # If 1, the timestep was found in the molanal.out file.
  '$found_volume',        # If 1, the volume was found in the molanal.out file.
  '$frame_count',         # The current frame read in.
  '$max_molecule_size',   # Molecules with more than this many atoms are a "polymer".
  '$min_mol_life',        # The minimum lifetime of a molecule in steps.
  '$poly_count',          # The number of polymer molecules found.
  '@poly_count_history',  # The history of the number of polymer molecules found.
  '$poly_numc',           # The number of carbons in the polymer.
  '$poly_numh',           # The number of hydrogens in the polymer.
  '$poly_numn',           # The number of nitrogens in the polymer.
  '$poly_numo',           # The number of oxygens in the polymer.
  '$poly_numf',           # The number of fluorines in the polymer.
  '@poly_numc_history',   # History variables for polymer properties.
  '@poly_numn_history',
  '@poly_numo_history',
  '@poly_numf_history',
  '@poly_numh_history',
  '$print_time_unit',     # The time unit conversion factor for printing histories.
  '$react_iter_max',      # Number of iterations to use in finding reactions.
  '$reaction_count_floor',# Only reactions with more than this many occurences have their
                          # Histories printed out.
  '$status',              # If 1, a frame was read successfully.
  '$stop_frame',          # Stop reading after this frame number.
  '$timestep',            # The time between frames in seconds.
  '%added_left',          # If 1, the indexed species was added to the left hand side of a reaction.
  '%added_right',         # If 1, the indexed species was added to the right hand side of a reaction.
  '%atoms',               # A hash of arrays containing atom labels for species.
  '%avgconc',             # The average concentration of each specie.
  '%charge',              # The average charge of a specie.
  '%conc_history',        # Per-timestep concentration history.
  '%dipole',              # Average dipole moment of a specie.
  '%found',               # If 1, the given molecule was found in the current frame.
  '%found_before',        # If 1, the given molecule was found in a previous frame
  '%life_count',          # Number of times a molecule has been created.
  '%lifetime',            # The average lifetime of a molecule.
  '%lifetime_histo',      # A histogram of lifetimes.
  '%mol_name',            # A hash to the name of a molecule, indexed by name + atom indices.
  '%molcount',            # The number of time a molecule has appeared in all frames.
  '%molecular_weight',    # The molecular weight of a molecule in the current frame, indexed by
                          # molecule name + tags.
  '$print_lifetime_histogram', # If non-zero, print out a histogram of species lifetimes.
  '$dq',                  # Spacing used in the charge histogram.
  '$qmax',                # Maximum charge used in the histogram.
  '%reaction_count',      # The number of times a reaction has occurred.
  '%reaction_history',    # The history of each reaction, giving how many times a reaction occurred as a
                          # function of time.
  '%reaction_lhs',        # The left hand side of a reaction.
  '%reaction_rhs',        # The right hand side of a reaction.
  '%rxn_flux',            # The reaction flux.
  '%total_atoms',         # The total number of atoms in a molecule, indexed by molecule name.
  '%total_charge',        # The average charge of a given molecule.
  '%total_charge_sq',     # The total q^2 of a given molecule.
  '%total_dipole',        # The average dipole moment of a given molecule.
  '%total_dipole_sq',     # The average dipole moment squared of a given molecule.
  '%total_molwt',         # The total molecular weight of a molecule, indexed by molecule name.
  '%total_time',          # The total time that a molecule was present.
  '@left',                # Molecules on the left hand side of a reaction.
  '@products',            # List of all molecules that have appeared on the current frame.
  '@reactants',           # List of all molecules that disappeared on the current frame.
  '@right',               # Molecules on the right hand side of a reaction.
  '@volume',               # The simulation cell volume on the current step, in angstroms^3.
  '@smp', 
  '$start_time' 
  ) ;
	     

sub read_frame  ;
sub add_left ;
sub add_right ;
sub check_atom_balance ;
sub make_reactions ;
sub print_totals ;
sub print_polymer_history ;
sub print_rates ;
sub print_conc_history ;
sub byflux ;
sub print_lifetime_histo ;

{
## This number is used to convert from the default input of volume in Ang.^3
## to moles/cc
    my $conc_fac = 1.0e24 / 6.02205e23 ;
    my $file ;
    my $i ;
    my $entry ;
    my %delete_entry ;
    my $name ;
    my $nqhisto ;
    my $life_histo_index ;
    my $atom_count ;
    my $old_atom_count ;

    $print_time_unit = 1.0e12 ;

    
    $frame_count = 0 ;
    if ( $#ARGV < 0 ) {
	print "Usage: findmolecules.pl <molanal.out>\n" ;
	exit(0) ;
    }
    $file = $ARGV[0] ;
    open(IN,"<$file") || die "Couldn't open $file\n" ;

    open(LOG,"> findmolecules.log") ;
    
    system ("date");
    print "\n";
    print "MOLECULE AND REACTION ANALYZER\n" ;
    print "By Larry Fried 11/06 \n\n" ;
    printf("Block averaging of %d steps\n", $block_size) ;
    printf("Molecules with more than %d atoms are treated as a polymer\n", $max_molecule_size) ;
    printf("Stop reading after frame %d (or end of file)\n", $stop_frame) ;
    printf("Molecules with concentrations less than %f will not be printed out\n", $conc_floor) ;
    printf("The minimum molecular lifetime is %f picoseconds\n", $min_mol_life) ;
    printf("The number of iterations to use in finding a reaction is %d\n", $react_iter_max) ;
    if ( $analyze_rates == 1 ) {
	printf("Reaction rates will be analyzed\n") ;
    } else {
	printf("Reaction rates will not be analyzed\n") ;
    }
    printf("Reaction histories will only be printed for reactions with more than %d occurrences\n",
	   $reaction_count_floor) ;

    $poly_count = 0.0 ;

##  Loop over frames in the molanal.out file.	

    do {

##      Reset the list of molecules in this frame.
	undef %found ;
	undef %charge ;
	undef %dipole ;
##	undef %mol_name ;

	if ( $frame_count > $stop_frame ) {
	    printf("Stopping at frame %d\n", $frame_count) ;
	    goto DONE ;
	} else {
	    print LOG "Reading frame $frame_count\n" ;
	}
##      Read in all molecules in this frame.
	read_frame() ;
	$block_count = int( ($frame_count-0.99999) / $block_size) ;
	$block_count1 = int( ($frame_count-1-0.99999) / $block_size) ;
	if ( $status != 1 ) {
	    printf("%d frames were read in\n", $frame_count) ;
	    goto DONE ;
	}
##      See which molecules are new, and which have been found before.
	@products = () ;
	@reactants = () ;
    	foreach $i ( keys (%found) ) {
	    if ( $found{$i} == 1 ) {
		## Update the lifetime of every molecule present.
		if ( $found_before{$i} == 1 ) {
		    $lifetime{$i}++ ;
		} else {
		    ## Start a new entry for molecules not found before.
		    $lifetime{$i} = 1 ;
		    $found_before{$i} = 1 ;
		    
		    $name = $mol_name{$i} ;
		    if ( ! defined $total_charge{$name} ) {
			$total_charge{$name} = 0.0 ;
			$total_charge_sq{$name} = 0.0 ;
			$total_dipole{$name} = 0.0 ;
			$total_dipole_sq{$name} = 0.0 ;
			$total_molwt{$name} = 0.0 ;
			$total_atoms{$name} = 0.0 ;
			$avgconc{$name} = 0.0 ;
		    }

		    ## Keep track of how many times this molecule has been "born".
		    $life_count{$mol_name{$i}}++ ;
		    
		    ## Add new molecules to the product list.
		    push(@products,$i) ; 
		}
		$molcount{$mol_name{$i}} += 1 ;

##              Normalize concentrations to the volume per step to support
##              variable volume simulations (NOT YET TESTED).	      
		$avgconc{$mol_name{$i}} += $conc_fac / $volume[$frame_count] ;

##
##              Accumulate averages of quantities of interest.
		$total_charge{$mol_name{$i}} += $charge{$i} ;
		
		$nqhisto = ($charge{$i} + $qmax)/$dq ;
		
		$charge_histo[$nqhisto] ++ ;
		$total_charge_sq{$mol_name{$i}} += $charge{$i} * $charge{$i} ;
		$total_dipole{$mol_name{$i}} += $dipole{$i} ;
		$total_dipole_sq{$mol_name{$i}} += $dipole{$i} * $dipole{$i} ;

		## Store the molecular weight by name.
                ## The weight could vary for lumped species, such as "polymer".
		$total_molwt{$mol_name{$i}} += $molecular_weight{$i} ;
		$total_atoms{$mol_name{$i}} += $#{ $atoms{$i} } + 1 ;
	    }
	}
	## Do a sanity check on the number of atoms.
	$atom_count = 0 ;
	foreach $i ( keys(%found) ) {
	    $atom_count += $#{ $atoms{$i} } + 1 ;
	}
	if ( $old_atom_count > 0 ) {
	    if ( $atom_count != $old_atom_count ) {
		die 'Failed atom count check\n' ;
	    }
	} else {
	    $old_atom_count = $atom_count ;
	}
##      Keep track of the concentration history of each variable
	foreach $i ( keys(%found) ) {
	    if ( $found{$i} == 1 ) {
#		print "CONC HISTORY: [$i] $conc_history{$mol_name{$i}}[$block_count] \n" ;
		$conc_history{$mol_name{$i}}[$block_count] += $conc_fac / $volume[$frame_count] ;
	    }
	}

##       Look for molecules that have disappeared
	foreach $i ( keys (%found_before) ) {
	    if ( $found_before{$i} == 1 ) {
		if ( $found{$i} == 0 ) {

#  		   The molecule was present before, but it isn't now.
#  		   Increment the total time this molecule was present.
		    $total_time{$mol_name{$i}} += $lifetime{$i} + $start_time ;
		    $life_histo_index = $lifetime{$i} + $start_time;
		    $lifetime_histo{$mol_name{$i}}[$life_histo_index] ++ ;

##                 Remember to free the memory associated with this molecule.
		    $delete_entry{$i} = 1 ;
		    push(@reactants,$i) ;
		} 
	    }
	}

	if ( $frame_count > 1 && $analyze_rates == 1 ) {
	    # Construct reactions by iterating between products and reactants, matching atoms,
	    # then balancing the equation.
	    make_reactions() ;

	}
##  Clean up unneeded memory.
	foreach $entry ( keys(%delete_entry) ) {
	    delete $found_before{$entry} ;
	    delete $mol_name{$entry} ;
	    delete $lifetime{$entry} ;
	    delete $atoms{$entry} ;
	    delete $molecular_weight{$entry} ;
	}
	undef %delete_entry ;
    } while ( $status == 1 ) ;
DONE: ;
#  Close out all found molecules.
    foreach $i ( keys (%found_before) ) {
	if ( $found_before{$i} == 1 ) {
	    $total_time{$mol_name{$i}} += $lifetime{$i} + $start_time ;
	}
    }

    ## Print totals of static properties.
    print_totals() ;

    ## Print rates averaged across the entire simulation.
    if ( $analyze_rates == 1 ) {
	print_rates() ;
    }
    
    ## Print out concentration histories.
    print_conc_history() ;

    ## Print out a histogram of lifetimes.
    if ( $print_lifetime_histogram == 1 ) {
	print_lifetime_histo() ;
    }

    ## Print out the history of the polymer composition.
    print_polymer_history() ;

    ## Print out reaction histories.
    if ( $analyze_rates == 1 ) {
	print_rate_history() ;
    }

}

sub read_frame {
## Sets variables:
#    $status: 1 on success, 0 otherwise
#    $number_molecules : Number of molecules found in this frame.
#    $found:  A list of all molecules found, indexed by tag.
#    $charge: A list of the charge of each molecule, indexed by tag.
#    $mol_name: A list of the name of each molecule, indexed by tag.
#    $molecular_weight: A list of the molecular weight of each molecule, indexed by tag.
#    $atoms: A list of all the atoms in the molecules found, indexed by tag.
#    $poly_count: The number of times a polymer species is found.
#    $found_timestep: Set to 1 if timestep found in the molanal output file.
#    $timestep: The time between frames read in.
    my ($numc, $numh, $numn, $numo, $numf) ;
    my @f ;
    my $found_number ;
    my $j ;
    my $index ;
    my $name ;
    my $thislist ;
    my $atoms_only ;
    my @namefields ;
    my $molwt ;
    my $ele ;
    my $i ;
    my @atom_list ;
    my $num ;
    my $elewt ;
    my $found_atoms ;
    my @tmp ;
    my $natoms ;
    my $q ;
    my $dip ;
    my $last_index ;
    my $found_start ;
    my $number_molecules ;

    $status = 1 ;
    $found_start = 0 ;
    while ( <IN> ) {
        if ( /Bonds require a lifetime of/ ) {
          @smp = split (" ") ;
          $start_time = $smp[5] ;
          print "min. lifetime (frames) = $start_time\n";
          $start_time = 0;
        }
	if ( /Beginning frame/ ) {
	    $found_start = 1 ;
	    last ;
	}
	if ( /The time between frames read in/ ) {
	    $found_timestep = 1 ;
	    @f = split(" ") ;
	    $timestep = $f[7] ;
	    if ( $timestep <= 0.0 || $timestep > 1.0e-10 ) {
		die "Bad timestep read in: $timestep\n" ;
	    }
	}
    }
    if ( $found_start == 0 ) {
	$status = 0 ;
	return ;
    } else {
	@f = split(" ", $_) ;
	$frame_count = $f[2] ;
	print LOG "Reading frame $frame_count\n" ;
    }
    $found_number = 0 ;
    $found_volume = 0 ;
    while ( <IN> ) {
	if ( /The box volume/ ) {
	    @f = split(" ") ;
	    $volume[$frame_count] = $f[4] ;
	    $found_volume = 1 ;
	}
	if ( /The number of molecules found/ ) {
	    $found_number = 1 ;
	    last ;
	}
    }
    if ( $found_number == 0 ) {
	$status = 0 ;
	return ;
    } else {
	@f = split(" ", $_) ;
	$number_molecules = $f[6] ;
    }
    if ( $found_volume == 0 ) {
	print "Warning: did not find the box volume: setting to 1.\n" ;
	$volume[$frame_count] = 1.0 ;
    }
    for ( $j = 0 ; $j < $number_molecules ; $j++ ) {
	$_ = <IN> ;
	if ( /Beginning molecule/ ) {
	    @f = split(" ", $_) ;
	    $index = $f[2] ;
	    if ( $index != $j + 1 ) {
		printf("Error: bad molecule index %d\n", $index) ;
		$status = 0 ;
		return ;
	    }
	}
	$_ = <IN> ;
	if ( /Name:/ ) {
	    chop $_ ;
	    $_ =~ s/Name: // ;
	    $name = $_ ;
	    $thislist = $name ;
	    $atoms_only = "" ;
	    undef @atom_list ;
#	    @atom_list = ("") ;
	    # Break the name up to parse the elemental composition.
	    @namefields = split(" ", $name) ;
	    $molwt = 0.0 ;
	    $numc = 0.0 ;
	    $numh = 0.0 ;
	    $numn = 0.0 ;
	    $numo = 0.0 ;
	    foreach $i ( @namefields ) {
#		print "Parsing field $i\n" ;
		$ele = substr($i,0,1) ;
		$num = substr($i,1) ;
#		print "Number of element $ele is $num\n" ;
		if ( $ele eq "C" ) {
		    $elewt = 12 ;
		    $numc = $num ;
		} elsif ( $ele eq "H" ) {
		    $elewt = 1 ;
		    $numh = $num ;
		} elsif ( $ele eq "F" ) {
		    $elewt = 19 ;
		    $numf = $num ;
		} elsif ( $ele eq "O" ) {
		    $elewt = 16 ;
		    $numo = $num ;
		} elsif ( $ele eq "N" ) {
		    $elewt = 14 ;
		    $numn = $num ;
		} elsif ( $ele eq "X" ) {
		    $elewt = 35.5 ;
		    $numn = $num ;
		} elsif ( $ele =~ /[0-9]/ ) {
		    last ;
		} else {
		    die "Bad element name $ele\n" ;
		}
		# Add up the molecule weight.
		$molwt += $elewt * $num ;
	    }
#	    print "The molecule weight is $molwt\n" ;
	} else {
	    $status = 0 ;
	    return ;
	}
	$found_atoms = 0 ;
	while ( <IN> ) {
	    if ( /Atom list/ ) {
		chop ;
		$_ =~ s/   Atom list:// ;
		$atoms_only = $atoms_only . $_ ;
		$thislist = $thislist . $_ ;
		$found_atoms = 1 ;
		@tmp = split(" ",$_) ;
#		print "ATOMS = ", @tmp, "\n" ;
		push(@atom_list,@tmp) ;
#		print "ATOM LIST = ", @atom_list, "\n" ;
	    } else {
		last ;
	    }
	}
	if ( $found_atoms == 0 ) {
	    $status = 0 ;
	    return ;
	}
	$natoms = $#atom_list ;
#	printf( "Size of the atom list = %d\n", $natoms ) ;
	if ( $natoms > $max_molecule_size ) {
	    $thislist = "polymer " . $atoms_only ;
	    $name = "polymer" ;
	    # Keep track of the polymer stoichiometry.
	    $poly_numc_history[$block_count] += $numc ;
	    $poly_numn_history[$block_count] += $numn ;
	    $poly_numo_history[$block_count] += $numo ;
	    $poly_numf_history[$block_count] += $numo ;
	    $poly_numh_history[$block_count] += $numh ;
	    $poly_count_history[$block_count] ++ ;
	    
	    $poly_numc += $numc ;
	    $poly_numn += $numn ;
	    $poly_numo += $numo ;
	    $poly_numf += $numf ;
	    $poly_numh += $numh ;
	    $poly_count++ ;
	} 
	$found{$thislist} = 1  ;
	if ( /Charge/ ) {
	    $found_charge = 1 ;
	    @f = split(" ") ;
	    $q = $f[2] ;
	    $charge{$thislist} = $q ;	
	    $_ = <IN> ;
	} 
	if ( /Dipole/ ) {
	    @f = split(" ") ;
	    $dip = $f[3] ;
	    $dipole{$thislist} = $dip ;	
	    $_ = <IN> ;
	} 
	if ( /Ending molecule/ ) {
	    @f = split(" ", $_) ;
	    $last_index = $f[2] ;
	    if ( $last_index != $index ) {
		printf("Bad index: %d\n", $last_index) ;
		$status = 0 ;
		return ;
	    }
## Use a tricky hash of arrays to hold all of the atoms for a molecule.
	    $atoms{$thislist} = [ @atom_list ] ;
#	    print "ATOMS($thislist) =", @{ $atoms{$thislist}}, " \n" ;
## See"Programming Perl", p. 268 Hashes of Arrays for syntax.	    
#	    foreach $ii ( 0 .. $#{ $atoms{$thislist} } ) {
#		print " atom $atoms{$thislist}[$ii] \n" ;
#	    }
	} else {
	    printf("Bad line: %s\n", $_) ;
	}
	$mol_name{$thislist} = $name ;
	$molecular_weight{$thislist} = $molwt ;
#	print "Molecule name $mol_name{$thislist} \n" ;
    }
    $status = 1 ;
    return ;
}

sub byconc 
# Used for sorting by concentration.
{
    $molcount{$b} <=> $molcount{$a}
}

sub bycount
# Used for sorting by reaction count
{
    $reaction_count{$b} <=> $reaction_count{$a}
}

sub byflux
# Used for sorting by reaction count
{
    $rxn_flux{$b} <=> $rxn_flux{$a}
}

sub add_right 
# Add product molecules to the right hand side of the reaction as needed.
{
##  Search across each reactant on the left hand side of the reaction.
    my $reac ;
    my $ridx ;
    my $ratm ;
    my $found_atom ;
    my $prod ;
    my $pidx ;
    my $patm ;
    
    foreach $reac ( @left ) {
	## Loop over each atom in the reactant.
	foreach $ridx ( 0 .. $#{ $atoms{$reac} } ) {
	    $ratm = $atoms{$reac}[$ridx] ;
	    $found_atom = 0 ;
	    
	    ## Search for the given reactant atom in each product.
	    foreach $prod ( @products ) {

		## Loop over each atom in the product.
		foreach $pidx ( 0 .. $#{ $atoms{$prod} } ) {
		    $patm = $atoms{$prod}[$pidx] ;

		    if ( $ratm == $patm ) {
			## Match found - between reactant and product atoms.

			if ( $added_right{$prod} != 1 ) {
			    ## If this product has not yet been added to the right hand side of
			    ## the reaction, add it.
			    push(@right,$prod) ;
			    $added_right{$prod} = 1 ;
			    $changed_right = 1 ;
#			    print "Added $prod to right\n" ;
			} 
			$found_atom = 1 ;
			goto FOUNDRATM ;
		    }
		}
	    }
	  FOUNDRATM:
#           Sanity check:  every reactant atom should be found in some product.
	    if ( $found_atom == 0 ) {
		die "Error: reactant and product list did not match\n" ;
	    } 
	}
    }
    return ;
}

sub add_left
# Add reactant molecules to the left side of the chemical reaction.
{
    my $reac ;
    my $ridx ;
    my $ratm ;
    my $found_atom ;
    my $prod ;
    my $pidx ;
    my $patm ;
    
    ## Search across each product on the right hand side of the reaction.
    foreach $prod ( @right ) {

	## Loop over each atom in the product.
	foreach $pidx ( 0 .. $#{ $atoms{$prod} } ) {
	    $patm = $atoms{$prod}[$pidx] ;
	    $found_atom = 0 ;

	    ## Search for the given product atom in each reactant.
	    foreach $reac ( @reactants ) {

		## Loop over each atom in the reactant.
		foreach $ridx ( 0 .. $#{ $atoms{$reac} } ) {
		    $ratm = $atoms{$reac}[$ridx] ;

		    if ( $ratm == $patm ) {
			## Match found - between reactant and product atoms.
			if ( $added_left{$reac} != 1 ) {
			    ## If this reactant has not yet been added to the left hand side of
			    ## the reaction, add it.
			    push(@left,$reac) ;
			    $added_left{$reac} = 1 ;
			    $changed_left = 1 ;
#			    print "Added $reac to left\n" ;
			} 
			$found_atom = 1 ;
			goto FOUNDPATM ;
		    }
		}
	    }
#  Sanity check:  every product atom should be found in some reactant.
	  FOUNDPATM:
	    if ( $found_atom == 0 ) {
		die "Error: reactant and product list did not match\n" ;
	    } 
	}
    }
}

sub check_atom_balance
# A last check to make sure all the atoms balance in the reaction.
{
    my $prod ;
    my $reac ;
    my $pidx ;
    my $ridx ;
    my $patm ;
    my $ratm ;
    my $found_atom ;


    ## Loop across each product molecule.
    foreach $prod ( @right ) {
	## Loop across every atom in the product.
	foreach $pidx ( 0 .. $#{ $atoms{$prod} } ) {
	    $patm = $atoms{$prod}[$pidx] ;
	    $found_atom = 0 ;
	    ## Each atom should be found in some reactant.
	    foreach $reac ( @left ) {
		foreach $ridx ( 0 .. $#{ $atoms{$reac} } ) {
		    $ratm = $atoms{$reac}[$ridx] ;
		    if ( $ratm == $patm ) {
			$found_atom = 1 ;
			goto FOUNDPATM ;
		    }
		}
	    }
	  FOUNDPATM:
	    if ( $found_atom == 0 ) {
		die "Error: reactant and product list did not match\n" ;
	    } 
	}
    }

    ## Loop across each reactant molecule.
    foreach $reac ( @left ) {
	## Loop across every atom in the reactant.
	foreach $ridx ( 0 .. $#{ $atoms{$reac} } ) {
	    $ratm = $atoms{$reac}[$ridx] ;
	    $found_atom = 0 ;
	    ## Each atom should be found in some product.
	    foreach $prod ( @right ) {
		foreach $pidx ( 0 .. $#{ $atoms{$prod} } ) {
		    $patm = $atoms{$prod}[$pidx] ;
		    if ( $ratm == $patm ) {
			$found_atom = 1 ;
			goto FOUNDRATM2 ;
		    }
		}
	    }
	  FOUNDRATM2:
#                 See if the desired reactant atom was found.
	    if ( $found_atom == 0 ) {
		die "Error: reactant and product list did not match\n" ;
	    } 
	}
    }
    return ;
}


sub make_reactions
## Analyze the given products and reactants to generate reactions.
{
    my $prod ;  # A product molecule.
    my $reac ; # A reactant molecule.
    my $iter ;  # The number of iterations.
    my $i ;     # Loop counter.
    my $l ;     # Entry from right.
    my $r ;     # Entry from left.
    my @rnames ;
    my @pnames ;
    my $reaction ;

    while ( $#products >= 0 && $products[0] ne "" ) {

	## Clean up variables for use below.
	undef @right ;
	undef @left ;
	undef %added_left ;
	undef %added_right ;

	$prod = $products[0] ;
	push(@right, $prod)  ;
	$added_right{$prod} = 1 ;

#               Start with a product and search for reactants with matching atoms.
#
	for ( $iter = 0 ; $iter <= $react_iter_max ; $iter++ ) {
	    $changed_right = 0 ;
	    $changed_left = 0 ;

	    add_left() ;
	    if ( $changed_left == 0 ) {
		last ;
	    }
	    add_right() ;
	    if ( $changed_right == 0 ) {
		last ;
	    }
	}
	if ( $iter <= $react_iter_max ) {
	    ## Make sure that the reaction satisfies detailed atom balance.
	    check_atom_balance() ;

##          Sort order of reactant and product names to avoid permutation problems.
	    @rnames = sort @left ;
	    @pnames = sort @right ;
	    undef $reaction ;
	    for ( $i = 0 ; $i <= $#rnames - 1 ; $i++ ) {
		$reaction = $reaction . $mol_name{$rnames[$i]} . " + " ;
	    }
	    $reaction = $reaction . $mol_name{$rnames[$#rnames]} . " => " ;

	    for ( $i = 0 ; $i <= $#pnames - 1 ; $i++ ) {
		$reaction = $reaction . $mol_name{$pnames[$i]} . " + " ;
	    }
	    $reaction = $reaction . $mol_name{$pnames[$#pnames]} ;
#	    print "REACTION: $reaction\n" ;
#	    print "RNAMES: @rnames \n" ;
#	    print "LEFT : @left\n" ;
#	    print "RIGHT: @right\n" ;
#	    exit(1) ;

	    print LOG "Reaction: $reaction\n" ;
	    $reaction_count{$reaction}++ ;
	    $reaction_history{$reaction}[$block_count1] += 1 ;
	    for ( $i = 0 ; $i <= $#rnames ; $i++ ) {
		$reaction_lhs{$reaction}[$i] = $mol_name{$rnames[$i]} ;
	    }
	    for ( $i = 0 ; $i <= $#pnames ; $i++ ) {
		$reaction_rhs{$reaction}[$i] = $mol_name{$pnames[$i]} ;
	    }

	} else {
	    die "COULD NOT FIND A CONSISTENT REACTION\n" ;
	}

	# Remove species from the right hand side of the reaction from the products list.
	foreach $l (@right ) {
	    for ( $i = 0 ; $i <= $#products ; $i++ ) {
		if ( $l eq $products[$i] ) {
		    splice(@products, $i, 1) ;
		}
	    }
	}

	# Remove species from the left hand side of the reaction from the reactants list.
	foreach $r (@left ) {
	    for ( $i = 0 ; $i <= $#reactants ; $i++ ) {
		if ( $r eq $reactants[$i] ) {
#			    print "Trimming $r from the reactants\n" ;
		    splice(@reactants, $i, 1) ;
		}
	    }
	}
    } 
    ## All products should have been put into some reaction.  Otherwise, we have an error.
    foreach $prod ( @products ) {
	die "ERROR REMAINING PRODUCT $prod\n" ;
    }

    ## All reactants should have been put into some reaction.  Otherwise, we have an error.
    foreach $reac ( @reactants ) {
	die "ERROR REMAINING REACTANT $reac\n" ;
    }

}

sub print_totals
## Print totals for static quantities.
{
    my $total_molecules ;
    my $i ;
    my $q2avg ;
    my $dipavg ;
    my $total_mass ;
    my $conc_sum ;
    my $mfrac_sum ;
    my $conc ;
    my $q ;
    my $dip ;
    my $molwt ;
    my $mfrac ;
    my $q_sq ;
    my $tst_sum ;
    my $tst_mfrac_sum ;
    my $time ;
    my $print_name ;
    my $delta_charge ;
    my $sum ;
    my $q ;
    my $j ;
    my $testsum ;
    my $dip2 ;
    my $dipfluc ;
    my $dip2avg ;

    $total_molecules = 0.0 ;
    foreach $i ( keys(%molcount) ) {
	$total_molecules += $molcount{$i} ;
    }
    $q2avg = 0.0 ;
    $dipavg = 0.0 ;

##  Calculate the total mass.
    $total_mass = 0.0 ;
    foreach $i ( sort byconc keys(%molcount) ) {
	## Mass = molar concentration * molecular weight.
  	$conc = $molcount{$i} ;
  	$conc = $conc / $total_molecules ;

	$molwt = $total_molwt{$i} / $molcount{$i} ;

  	$total_mass += $conc * $molwt ;
    }
##  Print out a table of molecular properties.	
    print("\nAverages for molecules\n") ;
    if ( $found_charge == 1 ) {
	print "                             Molecule type                            Mol. Frac. Mass Fr.  Mol. Wt. Avg. Q  Delta Q  Dipole  Delta Dip. Lifetime (ps)\n" ;
    } else {
	print "                             Molecule type                            Mol. Frac. Mass Fr.  Mol. Wt. Lifetime (ps) \n" ;
    }
    $mfrac_sum = 0.0 ;
    $conc_sum = 0.0 ;
    foreach $i ( sort byconc keys(%molcount) ) {
  	$conc = $molcount{$i} ;
  	$conc = $conc / $total_molecules ;
  	$q = $total_charge{$i} / $molcount{$i} ;
  	$dip = $total_dipole{$i} / $molcount{$i} ;
  	$dip2 = $total_dipole_sq{$i} / $molcount{$i} ;
	$dipfluc = sqrt($dip2 - $dip * $dip) ;
	$molwt = $total_molwt{$i} / $molcount{$i} ;
	$total_atoms{$i} /= $molcount{$i} ;

	$mfrac = $conc * $molwt / $total_mass ;
	$mfrac_sum += $mfrac ;
	$conc_sum += $conc ;
	$q_sq = $total_charge_sq{$i} / $molcount{$i} ;
##  	printf("Life count = $life_count{$i}\n") ;
  	if ( $life_count{$i} != 0 ) {
  	    $time = $total_time{$i} * $timestep * $print_time_unit /$life_count{$i};
#	    printf("Time = %f\n", $time) ;
  	} else {
  	    die "Error: zero life_count for $i\ n" ;
  	}
	$avgconc{$i} = $avgconc{$i} / $frame_count ;
	if ( $conc > $conc_floor && $time > $min_mol_life ) {
	    $print_name = substr($i,0,69) ;
	    if ( $found_charge == 1 ) {
		$delta_charge = sqrt($q_sq - $q * $q) ;
		printf("%70s   %7.4f %7.4f   %7.2f %7.2f %7.2f  %7.2f %7.2f    %7.4f\n", $print_name, $conc, $mfrac, $molwt, $q, $delta_charge, $dip, $dipfluc, $time) ;
	    } else {
		printf("%70s   %7.4f  %7.4f %7.2f %7.4f\n", $print_name, $conc, $mfrac, $molwt, $time) ;
	    }
	}
    }

##  Look for transition states and accumulate averages.
    if ( $found_timestep != 1 ) {
	printf("Warning: lifetime is reported in steps, not ps\n") ;
    }
    printf("\n\n") ;
    print("\nAverages for transition states\n") ;
    if ( $found_charge == 1 ) {
	print "                             Transition state                            Mol. Frac. Mass Fr.  Mol. Wt. Avg. Q  Delta Q  Dipole  Delta Dip. Lifetime (ps)\n" ;
    } else {
	print "                             Transition state                            Mol. Frac. Mass Fr.  Mol. Wt. Lifetime (ps)\n" ;
    }
    $mfrac_sum = 0.0 ;
    $conc_sum = 0.0 ;
    $q2avg = 0.0 ;
    $dipavg = 0.0 ;
    $dip2avg = 0.0 ;
    foreach $i ( sort byconc keys(%molcount) ) {
  	$conc = 1.0000001 * $molcount{$i} ;
  	$conc = $conc / $total_molecules ;
  	$q = $total_charge{$i} / $molcount{$i} ;
  	$dip = $total_dipole{$i} / $molcount{$i} ;
  	$dip2 = $total_dipole_sq{$i} / $molcount{$i} ;
	$dipfluc = sqrt($dip2 - $dip * $dip) ;
	$molwt = $total_molwt{$i} / $molcount{$i} ;
	$q_sq = $total_charge_sq{$i} / $molcount{$i} ;
	$mfrac = $conc * $molwt / $total_mass ;
	$mfrac_sum += $mfrac ;
	$conc_sum += $conc ;
##  	printf("Life count = $life_count{$i}\n") ;
  	if ( $life_count{$i} != 0 ) {
  	    $time = $total_time{$i} * $timestep * $print_time_unit / $life_count{$i};
#	    printf("Time = %f\n", $time) ;
  	} else {
  	    die "Error: zero life_count for $i\ n" ;
  	}

	if ( $time <= $min_mol_life ) {
#           Find concentration of transition states.
	    $tst_sum += $conc ;
	    $tst_mfrac_sum += $mfrac ;
	} else {
#           Accumulate average of molecular quantities.
	    $q2avg += $q_sq * $conc ;
	    $dipavg += $dip * $conc ;
	    $dip2avg += $dip2 * $conc ;
	}	    
	$print_name = substr($i,0,69) ;
	if ( $conc > $conc_floor && $time <= $min_mol_life ) {
	    if ( $found_charge == 1 ) {
		$delta_charge = sqrt($q_sq - $q * $q) ;
		printf("%70s   %7.4f %7.4f   %7.2f %7.2f %7.2f  %7.2f %7.2f    %7.4f\n", $print_name, $conc, $mfrac, $molwt, $q, $delta_charge, $dip, $dipfluc, $time) ;
	    } else {
		printf("%70s   %7.4f  %7.4f %7.2f %7.4f\n", $print_name, $conc, $mfrac, $molwt, $time) ;
	    }
	}
    }
    $q2avg /= (1.0 - $tst_sum) ;
    $dipavg /= (1.0 - $tst_sum) ;
    $dip2avg /= (1.0 - $tst_sum) ;

    if ( $poly_count > 0 ) {
	$poly_numc /= $poly_count ;
	$poly_numn /= $poly_count ;
	$poly_numo /= $poly_count ;
	$poly_numf /= $poly_count ;
	$poly_numh /= $poly_count ;
	printf("\nThe average polymer composition is : %5.2f H %5.2f C %5.2f N %5.2f O %5.2f F\n", $poly_numh, $poly_numc, $poly_numn, $poly_numo,
	       $poly_numf) ;
    }
    printf("\nThe mole fraction sum of transition states = %11.4f\n", $tst_sum) ;
    printf("The mass fraction sum of transition states = %11.4f\n", $tst_mfrac_sum) ;

    if ( $found_charge ) {
	printf("The <q^2>^1/2 for molecules = %11.4f e\n", sqrt($q2avg) ) ;
	printf("The <dipole> for molecules  = %11.4f D\n", $dipavg ) ;
	printf("The sqrt(<dipole^2>) for molecules  = %11.4f D\n", sqrt($dip2avg) ) ;
	printf("\nThe charge histogram\n") ;
	printf("q (e)   P(q)\n") ;
	
	$sum = 0.0 ;
	for ( $j = 0 ; $j <= $#charge_histo ; $j++ ) {
	    $sum += $charge_histo[$j] ;
	}
	
	$testsum = 0.0 ;
	for ( $j = 0 ; $j <= $#charge_histo ; $j++ ) {
	    $charge_histo[$j] /= $sum ;
	    $q = -$qmax + $j * $dq ;
	    $testsum += $q * $charge_histo[$j] ;
	    printf("%7.4f %7.4f\n", $q, $charge_histo[$j]) ;
	}
	if ( $testsum > 0.1 ) {
	    die 'Error: charge histogram sum test failed\n' ;
	}
    }
}

sub print_rates
## Print out reaction rates.
{
    my $i ;
    my $order ;
    my $ratecst ;
    my $k ;
    my $species ;
    my $aconc ;
    my $j ;
    my @reactions_to_print ;
    my @rhs ;
    my @lhs ;
    my $match ;
    my $l ;
    my %printed_rxn ;
    my $found_reverse ;
    my $rev ;
    my $printcount ;
    my $flux ;
    my $mass ;
    my %reverse_rxn ;

    print "\nReaction counts follow:\n" ;
    if ( $found_timestep == 1 ) {
	print "Rate constants in units of (moles/cc)^order / second\n" ;
    } else {
	print "Rate constants in units of (moles/cc)^order / frame\n" ;
    }
    print " Count    Rate Const. Reaction \n" ;
    foreach $i ( sort bycount keys(%reaction_count) ) {
	$order = $#{ $reaction_lhs{$i} } + 1 ;
	$ratecst = $reaction_count{$i} / $frame_count ;
	if ( $ratecst > 1.0e-12 ) {
	    for ( $k = 0 ; $k < $order ; $k++ ) {
		$species = $reaction_lhs{$i}[$k] ;
		$aconc = $avgconc{$species} ;
		if ( $aconc > 0.0 ) {
		    $ratecst /= (1.0000001 * $aconc) ;
		} else {
		    print "REACTION = $i\n" ;
		    print "ORDER = $order \n" ;
		    print "REACTION_COUNT = $reaction_count{$i} \n" ;
		    print "K = $k\n" ;
		    print "SPECIES = $species\n" ;
		    print "ACONC = $aconc \n" ;
		    print "RATECST = $ratecst \n" ;
		    die "Reactant species had zero concentration but rate not zero\n" ;
		}
	    }
	}
	if ( $found_timestep == 1 ) {
	    $ratecst = $ratecst / $timestep ;
	}
	printf("%5d   %11.4e  %s\n", $reaction_count{$i}, $ratecst, $i) ;
    }

    if ( $find_reverse_rxn == 1 ) {

	printf("\n\nReactions paired with reverse reactions\n") ;
	print " Count    Rate Const. Reaction \n" ;

	@reactions_to_print = sort bycount keys(%reaction_count)  ;
	$printcount = 0 ;

	foreach $i ( @reactions_to_print ) {
	    $order = $#{ $reaction_lhs{$i} } + 1 ;
	    $ratecst = $reaction_count{$i} / $frame_count ;
	    if ( $ratecst > 1.0e-12 ) {
		for ( $k = 0 ; $k < $order ; $k++ ) {
		    $species = $reaction_lhs{$i}[$k] ;
		    $aconc = $avgconc{$species} ;
		    if ( $aconc > 0.0 ) {
			$ratecst /= (1.0000001 * $aconc) ;
		    } else {
			print "REACTION = $i\n" ;
			print "ORDER = $order \n" ;
			print "REACTION_COUNT = $reaction_count{$i} \n" ;
			print "K = $k\n" ;
			print "SPECIES = $species\n" ;
			print "ACONC = $aconc \n" ;
			print "RATECST = $ratecst \n" ;
			die "Reactant species had zero concentration but rate not zero\n" ;
		    }
		}

	    }
	    if ( $found_timestep == 1 ) {
		$ratecst = $ratecst / $timestep ;
	    }
	    if ( $printed_rxn{$i} != 1 ) {
		$printcount++ ;
		printf("%d   %11.4e  %s\n", $reaction_count{$i}, $ratecst, $i) ;
		$printed_rxn{$i} = 1 ;
	    } else {
		next ;
	    }
	    $found_reverse = 0 ;
	    foreach $j ( @reactions_to_print ) {
		$match = 1 ;
		@lhs = @{$reaction_lhs{$i}} ;
		@rhs = @{$reaction_rhs{$j}} ;
		if ( $#lhs != $#rhs ) {
		    $match = 0 ;
		} else {
		    for ( $l = 0 ; $l <= $#lhs ; $l++ ) {
#		print "TEST $lhs[$l] $rhs[$l]\n" ;
			if ( $lhs[$l] ne $rhs[$l] ) {
			    $match = 0 ;
			    last ;
			}
		    }
		}
		if ( $match == 1 ) {
		    ## Explicitly match the reaction the other way around.
		    @lhs = @{$reaction_lhs{$j}} ;
		    @rhs = @{$reaction_rhs{$i}} ;
		    if ( $#lhs != $#rhs ) {
			$match = 0 ;
		    } else {
			for ( $l = 0 ; $l <= $#lhs ; $l++ ) {
#		print "TEST $lhs[$l] $rhs[$l]\n" ;
			    if ( $lhs[$l] ne $rhs[$l] ) {
				$match = 0 ;
				last ;
			    }
			}
		    }
		    if ( $match == 1 ) {
			$found_reverse = 1 ;
			$reverse_rxn{$i} = $j ;
			for ( $l = 0 ; $l <= $#lhs ; $l++ ) {
#		    print "MATCH $lhs[$l] $rhs[$l]\n" ;
			    ;
			}
		    }
		}
		if ( $match == 1 ) {
		    $rev = $j ;
		    $order = $#{ $reaction_lhs{$j} } + 1 ;
		    $ratecst = $reaction_count{$j} / $frame_count ;
		    if ( $ratecst > 1.0e-12 ) {
			for ( $k = 0 ; $k < $order ; $k++ ) {
			    $species = $reaction_lhs{$j}[$k] ;
			    $aconc = $avgconc{$species} ;
			    if ( $aconc > 0.0 ) {
				$ratecst /= (1.0000001 * $aconc) ;
			    } else {
				print "REACTION = $i\n" ;
				print "ORDER = $order \n" ;
				print "REACTION_COUNT = $reaction_count{$i} \n" ;
				print "K = $k\n" ;
				print "SPECIES = $species\n" ;
				print "ACONC = $aconc \n" ;
				print "RATECST = $ratecst \n" ;
				die "Reactant species had zero concentration but rate not zero\n" ;
			    }
			}
		    }
		    if ( $found_timestep == 1 ) {
			$ratecst = $ratecst / $timestep ;
		    }
		    if ( $printed_rxn{$j} != 1 ) {
			$printcount++ ;
			printf("%d   %11.4e  %s\n\n", $reaction_count{$j}, $ratecst, $j) ;
			$printed_rxn{$j} = 1 ;
		    }
		    last ;
		}
	    }
	    $order = $#{ $reaction_lhs{$i} } + 1 ;
	    $mass = 0 ;
	    for ( $k = 0 ; $k < $order ; $k++ ) {
		$species = $reaction_lhs{$i}[$k] ;
		$mass += $total_atoms{$species} ;
	    }
	    if ( $found_reverse == 1 ) {
		$flux = $mass * abs( $reaction_count{$rev} - $reaction_count{$i} ) ;
	    } else {
		$flux = $mass * $reaction_count{$i} ;
		print "\n" ;
#	    print "NO REVERSE FOUND\n" ;
	    }
	    $rxn_flux{$i} = $flux ;
	}

	printf("Flux       Reaction\n" ) ;
	foreach $i ( sort byflux keys(%rxn_flux) ) {

	    printf( "%7.2f   %s\n",$rxn_flux{$i}, $i) ;
	    if ( defined( $reverse_rxn{$i} ) ) {
		print "          ", $reverse_rxn{$i}, "\n" ;
	    }
	    print "\n" ;
	}

	if ( $printcount != $#reactions_to_print + 1 ) {
	    die "Error:  Did not print out the correct number of forward-reverse pairs\n" ;
	}
    }
}

sub print_polymer_history
# Print out the polymer concentration history
{
    my $i ;
    my $time ;

    printf("\nPolymer C History\n") ;
    for ( $i = 0 ; $i < $block_count ; $i++ ) {
	if ( $poly_count_history[$i] != 0 ) {
	    $poly_numc_history[$i] /= $poly_count_history[$i] ;
	}
	if ( $found_timestep == 1 ) {
	    $time = $i * $block_count * $timestep * $print_time_unit ;
	} else {
	    $time = $i * $block_count ;
	}
	printf("%7.3f %7.3f\n", $time, $poly_numc_history[$i] ) ;
    }

    printf("\nPolymer N History\n") ;
    for ( $i = 0 ; $i < $block_count ; $i++ ) {
	if ( $poly_count_history[$i] != 0 ) {
	    $poly_numn_history[$i] /= $poly_count_history[$i] ;
	}
	if ( $found_timestep == 1 ) {
	    $time = $i * $timestep * $print_time_unit ;
	} else {
	    $time = $i ;
	}
	printf("%7.3f %7.3f\n", $time, $poly_numn_history[$i] ) ;
    }

    printf("\nPolymer O History\n") ;
    for ( $i = 0 ; $i < $block_count ; $i++ ) {
	if ( $poly_count_history[$i] != 0 ) {
	    $poly_numo_history[$i] /= $poly_count_history[$i] ;
	}
	if ( $found_timestep == 1 ) {
	    $time = $i * $timestep * $print_time_unit ;
	} else {
	    $time = $i ;
	}
	printf("%7.3f %7.3f\n", $time, $poly_numo_history[$i] ) ;
    }


    printf("\nPolymer F History\n") ;
    for ( $i = 0 ; $i < $block_count ; $i++ ) {
	if ( $poly_count_history[$i] != 0 ) {
	    $poly_numf_history[$i] /= $poly_count_history[$i] ;
	}
	if ( $found_timestep == 1 ) {
	    $time = $i * $timestep * $print_time_unit ;
	} else {
	    $time = $i ;
	}
	printf("%7.3f %7.3f\n", $time, $poly_numf_history[$i] ) ;
    }

    printf("\nPolymer H History\n") ;
    for ( $i = 0 ; $i < $block_count ; $i++ ) {
	if ( $poly_count_history[$i] != 0 ) {
	    $poly_numh_history[$i] /= $poly_count_history[$i] ;
	}
	if ( $found_timestep == 1 ) {
	    $time = $i * $timestep * $print_time_unit ;
	} else {
	    $time = $i ;
	}
	printf("%7.3f %7.3f\n", $time, $poly_numh_history[$i] ) ;
    }
}

sub print_conc_history
## Print out the concentration history.
{
    my $i ;
    my $j ;
    my $k ;
    my $l ;
    my $m ;
    my $avg ;
    my $time ;

    print "\n" ;
    foreach $i ( sort byconc keys(%molcount) ) {
	if ( $avgconc{$i} > $conc_floor ) {
	    print "\nConcentration history for $i\n" ;
	    print "Time (ps)  Moles/cc\n" ;
	    for ( $j = 0 ; $j + $block_size < $frame_count ; $j+= $block_size ) {
		$m = $j / $block_size ;
		$conc_history{$i}[$m] = $conc_history{$i}[$m] / $block_size ;
		$avg = $conc_history{$i}[$m] ;
		if ( $found_timestep == 1 ) {
		    $time = $j * $timestep * $print_time_unit ;
		} else {
		    $time = $j ;
		}
		printf("%11.4e   %11.4e\n", $time, $avg ) ;
	    }
	}
    }
}


sub print_lifetime_histo
## Print out the lifetime histogram.
{
    my $i ;
    my $j ;
    my $k ;
    my $l ;
    my $m ;
    my $avg ;
    my $time ;
    my $n ;
    my $sum ;

    print "\n" ;
    foreach $i ( sort byconc keys(%molcount) ) {
	if ( $avgconc{$i} > $conc_floor ) {
	    print "\nLifetime histogram for $i\n" ;
	    print "Time (ps)  Count\n" ;
	    $n = $#{ $lifetime_histo{$i} } ;
	    $sum = 0.0 ;
	    for ( $j = 0 ; $j <= $n ; $j++ ) {
		$sum += $lifetime_histo{$i}[$j] ;
	    }
	    for ( $j = 0 ; $j <= $n ; $j++ ) {
		$avg = $lifetime_histo{$i}[$j] / $sum ;
		if ( $found_timestep == 1 ) {
		    $time = $j * $timestep * $print_time_unit ;
		} else {
		    $time = $j ;
		}
		printf("%11.4e   %11.4e\n", $time, $avg ) ;
	    }
	}
    }
}

sub print_rate_history
## Print out the reaction rate history.
{
    my $i ;
    my $j ;
    my $avg ;
    my $m ;
    my $k ;
    my $l ;
    my $order ;
    my $ratecst ;
    my $species ;
    my $aconc ;
    my $time ;

    print "\n" ;
    foreach $i ( sort bycount keys(%::reaction_count) ) {
	if ( $reaction_count{$i} > $reaction_count_floor ) {
	    print "\nReaction history for $i\n" ;
	    print "  Time        Avg. num. of rxns.      Rate constant\n" ;
	    $order = $#{ $reaction_lhs{$i} } + 1 ;
	    if ( $found_timestep == 1 ) {
		print "   (ps)          (1/ps)            1/(seconds * (mol/cc)^$order) \n" ;
	    } else {
		print "   (frames)      (1/frames)        1/(frames * (mol/cc)^$order) \n" ;
	    }
	    for ( $j = 0 ; $j + $block_size < $frame_count ; $j+= $block_size ) {
		$m = int($j / $block_size) ;
		$avg = $reaction_history{$i}[$m] / $block_size ;
		$ratecst = $avg ;
		if ( $ratecst > 1.0e-12 ) {
		    for ( $k = 0 ; $k < $order ; $k++ ) {
			$species = $reaction_lhs{$i}[$k] ;

			$aconc = $conc_history{$species}[$m] ;
			if ( $aconc > 0.0 ) {
			    $ratecst /= (1.0000001 * $aconc) ;
			} else {
			    die "Finite rate constant with zero concentration !\n" ;
			}
		    }
		}
		if ( $found_timestep == 1 ) {
		    $ratecst = $ratecst / $timestep ;
		    $time = $j * $timestep * $print_time_unit ;
		    $avg = $avg / ($timestep * $print_time_unit) ;
		} else {
		    $time = $j ;
		}
		printf("%11.4e       %11.4e           %11.4e\n", $time, $avg, $ratecst ) ;
	    }
	}
    }
    print "\n" ;

}

