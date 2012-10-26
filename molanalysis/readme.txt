Tools for the analysis of reactive molecular simulations (molanal)

Larry Fried  11/10/06

These are molecular analysis (molanal) tools.  They are used as follows:

1.  molanal.new - 
  
Usage:  molanal.new <trajectory file> <charge file> > molanal.out

This file takes as input the trajectory and charge (CHR.DAT) files
created by dftb_llnl.  The original DFTB does not accumulate charge
data.  The trajectory format is very simple, so it should
be easy to write a converter from another program.  The program will
loop over each frame of the trajectory, and report which molecules are
present in the simulation.  An xyz output file (molanal.xyz) is
created in which the molecules are "re-constructed" from the raw
trajectory.  If charges are given to the program, the dipole moment
wrt the center of mass and the net charge are reported for each
molecule.  molanal.new works for a non-orthorhombic unit cell.  Since
the code loops over a fixed number of periodic images, it could fail
to give correct results for a severely distorted cell.  For such
problems, the MAXWRAP parameter in the code should be increased.

bonds.dat:
Bond distances used to find molecules are taken from the file
bonds.dat.  As of 11/10/06, the bonds.dat file takes additional
input parameters.

First line: Minimum bond lifetime in MD steps.  Two atoms will not be
considered bonded unless they are continuously within the cutoff
distance for at least the given number of steps.  To get the behavior
of the previous code version, enter 0 here.  A minimum bond lifetime
of 10-30 fs seems reasonable.

New feature (NG, 11/07): In order to use different lifetimes for different 
bond types, set the minimum bond lifetime to be negative.  
Note: this feature requires the use of the file bond_times.dat (see below).

Second line: Number of copies for xyz file.  The specified number of
periodic copies of the simulation cell in each direction will be
printed to the xyz file.  This means that if you have natoms in the
cell, there will be natoms * (2*copies+1)^3 atoms output each step to
the xyz file.  Adding copies helps to make the trajectory look better.

Third line:  The time between frames in the gen file in seconds.

Fourth line:  If 1 is given, a directory called "molecules" will be
created containing a separate xyz file for each molecule found in the
simulation.  The file naming is similar to that used in the molanal.out
output file.  If 0 is given, no per-molecule xyz files will be generated.
The xyz files will reflect the last frame containing a particular molecule.

Additional lines:  Bond distances. Here is the format:
H H 1.1
H O 1.2

etc.

Here is a sample bonds.dat file:
-- Begin bonds.dat
20
1
2.41e-15
1
H H 0.900
C H 1.400
N H 1.450
O H 1.350
C C 1.900
C N 1.800
C O 1.800
N N 1.750
N O 1.650
O O 1.700
-- End bonds.dat

bond_times.dat: (optional)
Use this file to set different bond cut times for different atom pairs (see
above comment on bonds.dat).  The format is similar to the last part of 
bonds.dat. The bond lifetimes are given in terms of the MD time steps:
H H 15
H O 41

etc.

2.  gr_calc

Usage: gr_calc < gr_calc.in > gr_calc.out

This program calculates the radial distribution function for given element
pairs.  The program is written in fortran 90, so you will need an f90 compiler.
g95 or gfortran are freely available.

Here is an example gr_calc.in file:

-- Begin gr_calc.in
tatb_113cj5.out  ! Trajectory file.
4        ! input format type (4 = DFTB)
4 4      ! Desired atom types.
0.05     ! dr
2        ! number of super-cells added in each direction. (use 2 for TATB).
100.0    ! time blocking	
1        ! number of time blocks.
0.00025  ! time step.
10       ! Diffusion time block.
-- End gr_calc.in

The first line is the trajectory file.  The second line is the input
format type.  ONLY DFTB (4) HAS BEEN TESTED.  The third line gives the
desired atom types.  So if the DFTB trajectory file lists

"C" "N" "O" "H"

Then 1 is "C", 2 is "N", 3 is "O" and 4 is "H". The example file
calculates the H-H radial distribution file.  The dr parameter is
the bin size for the radial distribution function.  The gr_calc program
works with non-orthorhombic unit cells.  For such cells, a super-cell
search is done to determine distances.  Convergence with respect to 
the number of super-cells should be checked.  1 is adequate for 
orthorhombic cells.  

The time blocking parameter gives the option of calculating G(r) for
multiple time slices of the trajectory.  Give a large number if you
just want one time slice.  The number of time blocks desired is then
given.  The time step parameter is the time between configurations in
the trajectory file.  This is used in calculating the time blocks.

The diffusion time block is used in calculating the diffusion
coefficient and mean squared displacement.  The displacement of each
particle is placed into a bin with a width given (in steps) by the
diffusion time block.  Also, every nth time step is used for a new
initial condition in calculating the displacement, where n is also the
diffusion time block.  Diffusion constants are calculated by the
formula D = 1/6 d/dt < |r(t) - r(t0)| >^2.  The trajectory is divided
into 5 segments.  Each segment is used to calculate the time
derivative in the formula through finite differences.  This leads to
some scatter in the estimates, which is realistic given the difficulty
in calculating accurate diffusion constants.  The 1st estimate of the
diffusion constant starts from 0 time and thus includes the "inertial"
part of the system response.  It should not be used in reported
diffusion constants.

One reason for calculating the g(r) is to estimate a bond distance.
Two estimates are automatically printed out:  the second value where
g(r) = 1, and the first minimum of g(r) where g(r) < 1.  The second
was used by Goldman et al. in studying water.  The first was suggested
in a recent PRL on high temperature water.  Note that in some cases
one or both of the criteria can fail to give an answer.  In that case
there may not be a covalent bond between the species involved, or the
fraction of bonds may be very low.

3.  findmolecules.pl 

Usage: findmolecules.pl molanal.out > findmolecules.out

This program takes the output of molanal.new, and reports the 
average concentration, lifetime, dipole moment, and charge of each species.
The program has been modified as of 11/10/06 to calculate more quantities.
The code now also defines a "polymer" composite species.  This is useful
for simulations where a great number of polymeric species are formed.
All these species are now lumped into a single "polymer".  There
are a number of user-defined inputs inside the findmolecules.pl program.
These inputs need to be edited before running the program.

Here is the section of the program to edit:
##################################
#  USER-DEFINED VARIABLES
#
# The block_size is the number of steps to average concentrations over.
$block_size = 100 ;
# Molecules over this size will be counted as a "polymer".
$max_molecule_size = 10 ;
# Stop reading after this many frames. Or when end of file is found.
$stop_frame = 100000 ;
# Do not print out molecules with concentration less than this (in units of mol/cc).
$conc_floor = 1.0e-04 ;
# Minimum lifetime for a molecule (in ps).
$min_mol_life = 0.02 ;
# Number of iterations to use in finding reactions.
$react_iter_max = 10 ;
# Set to 1 in order to analyze rates, 0 otherwise.
$analyze_rates = 1 ;
# Reaction histories will be printed only for reactions with more than 
# this number of reaction events.
$reaction_count_floor = 10 ;
# If set to 1, calculate reverse reaction pairs.  This is somewhat slow right now.
$find_reverse_rxn = 0 ;
#  END USER-DEFINED VARIABLES
##################################

The block_size determines how many steps are used in averaging output
concentration histories.  

Some tuning of max_molecule_size is desirable.  Any molecule with more
then max_molecule_size is counted as a "polymer".  This prevents an
explosion of complexity in simulations where very many polymeric species
are formed.  If max_molecule_size is set to a large number (greater than
the number of atoms in the simulation), no "polymer" species will be created
and all molecules will be tracked.

The trajectory file will only be read up to the given stop_frame.  This is
useful for short exploratory calculations while determining parameters.

Only molecules with concentration larger than conc_floor are printed out.

A bonded entity is counted as a "molecule" if is lasts longer than the
number of steps given by min_mol_life.  Otherwise, it is counted as a
"transition state".

If analyze_rates is set to 1, reaction rates will be found as well as
reaction histories.  If it is set to 0, no reaction rate analysis will be
performed.

findmolecules.pl prints the each frame it reads into findmolecules.log.
It is somewhat helpful to look at this file if the program takes a long
time to run.  The program run time is sensitive to the max_molecule_size
for some simulations.

findmolecules.pl has been updated as of 11/10/06 to print out 
reactions as well as species.  Reactions are identified based on
the loss or gain of molecules between two frames.  A new file
findmolecules.log, is created when the program is run.  This
file has a list of which reactions occur during a particular frame.

If find_reverse_rxn is set to 1, reverse reactions are found and the net
reactive flux for the forward and reverse reactions are determined.  The 
reactive flux for a reaction i is defined as natoms * |nforward - nreverse|,
where natoms is the number of atoms on one side of the reaction.  nforward
is the number of forward reactions, and nreverse is the number of reverse
reactions.

4.  molanal.pl  

Usage: molanal.pl <gen name> <chr name>. 

This program runs molanal.new on a set of numbered trajectory and
charge files.  Since dftb_llnl no longer creates numbered files, this
program is mostly obsolete.  THIS HAS NOT BEEN TESTED IN A LONG TIME.

5.  gentoarc.pl

Usage: gentoarc.pl < <gen file>

This program reads in dftb trajectory (gen) file, and creates a
corresponding MSI arc file.  This file can be read in by Materials
Studio or by Cerius2.


6. molanal

This is the "classic" version of the molecular analyzer.  It only
works for orthorhombic unit cells.  It does not analyze charges.  Bond
lengths are compiled in.  Currently it is faster than molanal.new.
THIS HAS NOT BEEN TESTED IN A LONG TIME.  molanal.new supports bond
lifetimes, which gives a much better bond definition.

7. cpmdtogen.pl

This program converts a CPMD trajectory to DFTB gen format. The 
trajectory must be contiguous (no overlapping restart data).
The gen format can then be read by other tools in this directory.

Usage:  cpmdtogen.pl <TRAJECTORY> <cpmddat> <trajectory.gen>

The TRAJECTORY file is the cpmd trajectory dump.  The cpmddat
file contains extra data used to create the gen file.  Here
is an example:

- begin cpmd.dat -
7.9273 0.0 0.0
0.0 7.9273 0.0
0.0 0.0 7.9273
F 54
H 54
- end cpmd.dat -

first three lines contain the unit cell vectors.  They do not have to be
orthogonal.  The next lines contain the element types and counts used
in the trajectory file.  In the example, there are 54 HF molecules.  The
first 54 atoms are fluorines and the next 54 atoms are hydrogens.

The file trajectory.gen is the DFTB gen-format file created as output.

