/********************************************************************

A simple molecular analyzer for DFTB simulations.

usage: molanal.new <gen file> <chr file>

The chr (charge) file is optional.
The gen (coordinate) file is mandatory.

Non-orthorhombic simulation cells are supported.

Output:  Molecules found are printed out with total charges,
         if supplied.
         An xyz file suitable for jmol is created in the
	 file molanal.xyz.

Bond distances are now read from the file bonds.dat.

Larry Fried, 5/30/2005.

A bond lifetime criterion is now implemented.  Fixed a bug in 
dipole calculations when wrap_com was not called.

Larry Fried, 11/10/06

*********************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

#define BUFSIZE 1024  /** Size of buffer for IO **/
#define MAXATOM 820  /** Maximum number of atoms **/
#define MAXWRAP 2   /* Maximum number of lattice vectors to search over.
		     Use a larger number for a highly non-orthogonal cell.
		    */
#define MAXELE 6    /** Maximum number of elements **/
#define MAXHIST 150  /** Maximum number of history steps for bond determination. */
#define NOTSET 9999   /* Number to identify whether history variables have been set. */

int not_in_molecule(int isearch , int mol_list[MAXATOM][MAXATOM],
		    int nmolecule) ;
void find_molecules(double *x, double *y, double *z, 
  		    double a[3], double b[3], double c[3],
		    int mol_list[MAXATOM][MAXATOM], 
		    int natom, int *nmolecule,
		    int *type, char **element,
		    int bond_list[MAXATOM][MAXATOM],
                    int bond_time[MAXELE][MAXELE] ) ;
int is_bonded(int j, int k, double *x, double *y, double *z,
	      double a[3], double b[3], double c[3],
	      int *type, char **element, int bond_time[MAXELE][MAXELE]) ;
void find_molecule(int j, 
		   double *x, double *y, double *z,
		   double a[3],  double b[3], double c[3],
		   int mol_list[MAXATOM][MAXATOM],
		   int natom,
		   int nmolecule,
		   int *type,
		   char **element,
		   int bond_list[MAXATOM][MAXATOM],
                   int bond_time[MAXELE][MAXELE] ) ;
void wrap_atoms(double *x, double *y, double *z, 
		double a[3], double b[3], double c[3],
		int bond_list[MAXATOM][MAXATOM], 
		int natom)  ;
double wrap_atom(int j, int k, double *x, double *y, double *z,
		 double a[3], double b[3], double c[3], int maxwrap,
		 int do_wrap) ;
void wrap_molecules(double *x, double *y, double *z, 
		    double a[3], double b[3], double c[3],
		    int mol_list[MAXATOM][MAXATOM], int nmolecule,
		    int bond_list[MAXATOM][MAXATOM], int natom ) ;
void wrap_com(double *x, double *y, double *z, 
	      double a[3], double b[3], double c[3],
	      int mol_list[MAXATOM][MAXATOM], int nmolecule, int natom,
	      double rcomx[MAXATOM], double rcomy[MAXATOM], 
	      double rcomz[MAXATOM], int type[MAXATOM],
	      char *element[MAXELE]);
void center_all(double *x, double *y, double *z, 
	      double a[3], double b[3], double c[3],
	      int mol_list[MAXATOM][MAXATOM], int nmolecule) ;
double r2bond(char *ea, char *eb) ;
int nint(double num) ;
void print_molecule(int mol_list[MAXATOM][MAXATOM],
		    int bond_list[MAXATOM][MAXATOM],
		    int *type, char **element, int mol, 
		    int read_charge,
		    double q[MAXATOM],
		    double x[MAXATOM], double y[MAXATOM], double z[MAXATOM],
		    double rcomx[MAXATOM], double rcomy[MAXATOM],
		    double rcomz[MAXATOM], int printed_atom[MAXATOM],
		    int dump_structures)  ;
void inversebox(double a[3], double b[3], double c[3],
		double invbox[3][3])  ;
static void wrap_in_box(double a[3], double b[3], double c[3],
			double invbox[3][3],
			int natom, double x[MAXATOM], double y[MAXATOM], 
			double z[MAXATOM])  ;
void wrap_pairs(double *x, double *y, double *z, 
		double a[3], double b[3], double c[3],
		int mol_list[MAXATOM][MAXATOM], int nmolecule,
		int bond_list[MAXATOM][MAXATOM], int natom, int j, int k,
		int wrapped[MAXATOM]) ;
int ele_index(char *ea);
double atm_mass(char *ea) ;
void read_bonds(double rbond[MAXELE][MAXELE], int *duration, int *xyz_copies, double *time_step,
		int *dump_structures, int bond_time[MAXELE][MAXELE]) ;
double boxvol(double a[3], double b[3], double c[3]) ;


double xold[MAXHIST][MAXATOM], yold[MAXHIST][MAXATOM], zold[MAXHIST][MAXATOM] ;
int bond_duration = 21 ; /* Determines how many frames a bond must exist. */
double rbond[MAXELE][MAXELE] ;
int bond_time[MAXELE][MAXELE] ;

int main(int argc, char** argv) 
{
  char buf[BUFSIZE], buf1[BUFSIZE][BUFSIZE] ;
  double x[MAXATOM], y[MAXATOM], z[MAXATOM], q[MAXATOM] ;
  double rcomx[MAXATOM], rcomy[MAXATOM], rcomz[MAXATOM] ;
  double a[3], b[3], c[3] ;
  double invbox[3][3] ;
  int nmolecule, do_wrap_atoms ;
  int natom, i, j, k, type[MAXATOM], mol_list[MAXATOM][MAXATOM],
    bond_list[MAXATOM][MAXATOM] ;
  char *element[MAXELE], tmp_ele[MAXELE][MAXELE] ;
  FILE *fout ;
  FILE *fgen ;    /* The gen (coordinate) file. */
  FILE *fchr ;    /* The chr (charge) file */
  int read_charge = 0 ; /* Whether charges should be read. */
  int natom1 ;    /* Number of atoms in the chr file */
  double dx, dy, dz ;
  int xyz_copies = 0 ;  /* Print this many copies in each direction for the xyz file. */
  int frame = 1 ;
  int ix, iy, iz ;
  int printed_atom[MAXATOM] ;  /* Whether an atom was printed out. */
  double vbox ;        /* The volume of the box. */
  double time_step ;    /* The time step between saved frames. */
  int dump_structures = 0 ;  /* Whether to dump structures. */
  int bond_time_max ;

  do_wrap_atoms = 1 ;  /** Whether to wrap atoms. This is required for correct
			molecule identification. **/

  read_bonds(rbond, &bond_duration, &xyz_copies, &time_step, &dump_structures, bond_time) ;

  if ( dump_structures ) {
    system("rm -rf molecules") ;
    system("mkdir molecules") ;
  }

  bond_time_max = -1 ;
  for ( i = 0 ; i < MAXELE ; i++ ) {
    for ( j = i + 1 ; j < MAXELE ; j++ ) {
      rbond[i][j] = rbond[j][i] ;
      bond_time[i][j] = bond_time[j][i] ;
      if (bond_time[i][j] > bond_time_max) {
        bond_time_max = bond_time[i][j] ;
      }
    }
  }
  printf ("max. bond time cutoff (time steps) = %d\n", bond_time_max);

  if ( bond_time_max >= MAXHIST ) {
    printf("Error: bond duration is too big.  Recompile with a bigger MAXHIST\n") ;
    exit(1) ;
  }

  for ( i = 0 ; i < MAXATOM ; i++ ) {
    q[i] = 0.0 ;
    for ( j = 0 ; j < bond_time_max ; j++ ) {
      xold[j][i] = NOTSET ;
      yold[j][i] = NOTSET ;
      zold[j][i] = NOTSET ;
    }
  }
  if ( argc >= 2 ) {
    printf ("Trajectory file: %s\n",argv[1]);
    fgen = fopen(argv[1],"r") ;
    if ( fgen == NULL ) {
      printf("Error: could not open file %s\n", argv[1]) ;
      exit(1) ;
    }
  } else {
    printf("Error: The input file must be specified\n") ;
    exit(1) ;
  }
  if ( argc == 3 ) {
    read_charge = 1 ;
    printf ("Charge file: %s\n",argv[2]);
    fchr = fopen(argv[2],"r") ;
    if ( fchr == NULL ) {
      printf("Error: could not open file %s\n", argv[2]) ;
      exit(1) ;
    } 
  }
  fout = fopen("molanal.xyz", "w") ;

  while ( ! feof(fgen) || ferror(fgen) ) {
    /** Defensive programming. */
    for ( j = 0 ; j < MAXATOM ; j++ ) {
      x[j] = 0.0 ;
      y[j] = 0.0 ;
      z[j] = 0.0 ;
      q[j] = 0.0 ;
      rcomx[j] = -100.0 ;
      rcomy[0] = -100.0 ;
      rcomz[j] = -100.0 ;
      type[j] = 0 ;
    }
    for ( j = 0 ; j < MAXELE ; j++ ) {
      element[j] = '\0' ;
    }
    if ( argc == 3 ) {
      for ( i = 0 ; i < 4 ; i++ ) {
	fgets(buf,BUFSIZ,fchr) ;
	if ( feof(fchr) ) {
	  printf("End of charge file reached\n") ;
	  exit(0) ;
	}
	if ( buf[0] != '#' ) {
	  printf("Bad format in charge file.\n") ;
	  exit(1) ;
	}
      }
      fscanf(fchr,"#   %d\n", &natom1) ;
      for ( i = 0 ; i < natom1 ; i++ ) {
	fgets(buf,BUFSIZ,fchr) ;
	if ( sscanf(buf,"%d  %lf", &j, &q[i]) != 2 ) {
	  printf("Error: could not read in a charge\n") ;
	  exit(1) ;
	} 
      }
      if ( feof(fchr) || ferror(fchr) ) {
	printf("Error: end of CHR file found\n") ;
	exit(1) ;
      }
    }
    while ( fgets(buf,BUFSIZ,fgen) ) {
      if ( buf[0] != '#' ) {
	break ;
      }
    }
    if ( feof(fgen) || ferror(fgen) ) {
      break ;
    }
    sscanf(buf,"%d", &natom) ;
    if ( read_charge && natom != natom1 ) {
      printf("Error: number of atoms in the chr and gen files are different\n") ;
      exit(1) ;
    }
    if ( natom > MAXATOM ) {
      printf("Error: the number of atoms is too large.  Increase MAXATOM and recompile\n") ;
    }
    fgets(buf,BUFSIZ,fgen) ; 
    for ( j = 0 ; buf[j] != 0 ; j++ ) {
      if ( buf[j] == '\'' ) {
	buf[j] = ' ' ;
      }
    }
    /** Parse the input line, get the element types **/
    k = 0 ;
    for ( j = 0 ; buf[j] ; j++ ) {
      if ( isalpha((int)buf[j]) ) {
        if ( isalpha((int)buf[j+1]) ) {
          tmp_ele[k][0] = buf[j];
          tmp_ele[k][1] = buf[j+1];
          tmp_ele[k][2] = '\0';
          element[k] = tmp_ele[k];
          j++;
        } else {
          tmp_ele[k][0] = buf[j];
          tmp_ele[k][1] = '\0';
          element[k] = tmp_ele[k];
        }
        k++ ;
      }
    }
    for ( j = 0 ; j < natom ; j++ ) {
      fgets(buf,BUFSIZ,fgen) ;
      sscanf(buf,"%d %d %lf %lf %lf", &i, &type[j], &x[j], &y[j], &z[j]) ;
    }
    fgets(buf,BUFSIZ,fgen) ;
    
    fgets(buf,BUFSIZ,fgen) ;
    sscanf(buf,"%lf %lf %lf", &a[0], &a[1], &a[2]) ;
    fgets(buf,BUFSIZ,fgen) ;
    sscanf(buf,"%lf %lf %lf", &b[0], &b[1], &b[2]) ;
    fgets(buf,BUFSIZ,fgen) ;
    sscanf(buf,"%lf %lf %lf", &c[0], &c[1], &c[2]) ;
    /*    fgets(buf,BUFSIZ,fgen) ; */
    if ( a[0] <= 0.0 || b[1] <= 0.0 || c[2] <= 0.0 ) {
      printf("Bad cell size read in\n") ;
      exit(1) ;
    }
    /** Wrap all coordinates back into the primitive unit cell. */
    inversebox(a,b,c,invbox) ;
    wrap_in_box(a,b,c,invbox,natom,x,y,z) ;
    vbox = boxvol(a,b,c) ;

    for ( j = 0 ; j < MAXATOM ; j++ ) {
      for ( i = 0 ; i < MAXATOM ; i++ ) {
	mol_list[i][j] = - 1 ;
	bond_list[i][j] = 0 ;
      }
    }
    nmolecule = 0 ;
    /** Each atom is assigned to a molecule **/
    find_molecules(x,y,z,a,b,c,mol_list,natom,&nmolecule,
		   type,element,bond_list,bond_time) ;

    if ( frame > bond_time_max ) {
      printf("Beginning frame %d\n", frame - bond_time_max ) ;
      printf("The box volume = %11.7e\n", vbox) ;
      printf("The number of molecules found = %d\n", nmolecule) ;
    }
    
    /** Wrap each atom to make whole molecules **/
    /** wrap_atoms(x, y, z, a, b, c, bond_list, natom) ; **/
    wrap_molecules(x, y, z, a, b, c, mol_list, nmolecule,bond_list, natom) ;
    
    /** Wrap the center of mass of each molecule to be inside the box
     **/
    wrap_com(x,y,z,a,b,c,mol_list,nmolecule, natom,
	     rcomx, rcomy, rcomz,type,element);  

    if ( frame > bond_time_max ) {
      for ( j = 0 ; j < natom ; j++ ) {
	printed_atom[j] = 0 ;
      }
      for ( j = 0 ; j < nmolecule  ; j++ ) {
	printf("Beginning molecule %d\n",j+1) ;
	print_molecule(mol_list,bond_list,type,element,j,read_charge,q,
		       x,y,z,rcomx,rcomy,rcomz,printed_atom,dump_structures) ;
	printf("Ending molecule %d\n",j+1) ;
      }
      for ( j = 0 ; j < natom ; j++ ) {
	if ( printed_atom[j] == 0 ) {
	  printf("Error: did not print out atom %d\n", j) ;
	  exit(1) ;
	}
      }
      fprintf(fout,"%d\n\n", natom * (2*xyz_copies+1) * (2*xyz_copies+1) * (2*xyz_copies+1)) ;
      for ( ix = -xyz_copies; ix <= xyz_copies; ix++ ) {
	for ( iy = -xyz_copies; iy <= xyz_copies; iy++ ) {
	  for ( iz = -xyz_copies; iz <= xyz_copies; iz++ ) {
	    dx = ix * a[0] + iy * b[0] + iz * c[0] ;
	    dy = ix * a[1] + iy * b[1] + iz * c[1] ;
	    dz = ix * a[2] + iy * b[2] + iz * c[2] ;
	    for ( j = 0 ; j < natom ; j++ ) {
	      fprintf(fout,"%s %f %f %f\n", element[type[j]-1], 
		      dx + x[j], dy + y[j], dz + z[j]) ;
	    }
	  }
	}
      }
      printf("Ending frame %d\n", frame - bond_time_max) ;
    }

    if ( feof(fgen) || ferror(fgen) ) {
      break ;
    }
    for ( j = 0 ; j < natom ; j++ ) {
      for ( i = 0 ; i < bond_time_max - 1 ; i++ ) {
	xold[i][j] = xold[i+1][j] ;
	yold[i][j] = yold[i+1][j] ;
	zold[i][j] = zold[i+1][j] ;
      }
      if ( bond_time_max > 0 ) {
	xold[bond_time_max-1][j] = x[j] ;
	yold[bond_time_max-1][j] = y[j] ;
	zold[bond_time_max-1][j] = z[j] ;
      }
      /** 
      if ( j == 0 ) {
	for ( i = 0 ; i < bond_duration ; i++ ) {
	  printf("Xold[%d][0] = %f\n", i, xold[i][0]) ;
	}
      }
      **/
    }
    frame++ ;

  }
  fclose(fout) ;
  return(0) ;
}


void find_molecules(double *x, double *y, double *z, 
  		    double a[3], double b[3], double c[3],
		    int mol_list[MAXATOM][MAXATOM], 
		    int natom, int *nmolecule,
		    int *type, char **element, 
                    int bond_list[MAXATOM][MAXATOM],
                    int bond_time[MAXELE][MAXELE]) 
{
  int j, k, l ;
  int sort_list[MAXATOM] ;
  int used[MAXATOM] ;
  int best_idx ;
  int best ;

  for ( j = 0 ; j < natom ; j++ ) {
    if ( not_in_molecule(j,mol_list,*nmolecule) ) {
      mol_list[*nmolecule][0] = j ; 
      find_molecule(j,x,y,z,a,b,c,mol_list,natom,*nmolecule,
		    type, element,bond_list,bond_time) ;
#if(0)
 {  
   int k ;
      for ( k = 0 ; mol_list[*nmolecule][k] >= 0 ; k++ ) {
	printf("mol:%d %d %d\n", *nmolecule, k, mol_list[*nmolecule][k]) ;
      }
 }
#endif
      (*nmolecule)++ ;
    }
  }
  /* Simple sort so that the atom list is unique. */
  for ( j = 0 ; j < *nmolecule ; j++ ) {
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      used[k] = 0 ;
    }
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      best = 999999 ;
      best_idx = 0 ;
      for ( l = 0 ; mol_list[j][l] >= 0 ; l++ ) {
	if ( mol_list[j][l] < best &&
	     used[l] == 0 ) {
	  best = mol_list[j][l] ;
	  best_idx = l ;
	}
      }
      sort_list[k] = mol_list[j][best_idx] ;
      used[best_idx] = 1 ;
    }
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      mol_list[j][k] = sort_list[k] ;
    }
  }
}



void wrap_atoms(double *x, double *y, double *z, 
		double a[3], double b[3], double c[3],
		int bond_list[MAXATOM][MAXATOM], 
		int natom) 
     /* Wrap bonded atoms so that they are as close as possible to each other. */
{
  int j, k ;
  for ( j = 0 ; j < natom ; j++ ) {
    for ( k = j+1 ; k < natom ; k++ ) {
      if ( bond_list[j][k] == 1 ) {
	/**
	printf("Wrapping atoms %d and %d, x = %f, y = %f, z = %f\n",
	       j, k, x[k], y[k], z[k]) ;
	**/
	wrap_atom(j, k, x,y,z,a,b,c,MAXWRAP,1) ;
	/**
	printf("Done wrapping atoms %d and %d, x = %f, y = %f, z = %f\n",
	       j, k, x[k], y[k], z[k]) ;
	**/
      }
    }
  }
}


int not_in_molecule(int isearch , int mol_list[MAXATOM][MAXATOM],
		    int nmolecule)
     /** Return 1 if atom isearch is not in the given molecule, 0
	 otherwise **/
{
  int j, k ;
  for ( j = 0 ; j < nmolecule ; j++ ) {
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      if ( mol_list[j][k] == isearch ) {
	return ( 0 ) ;
      }
    }
  }
  return (1) ;
}

void find_molecule(int j, 
		   double *x, double *y, double *z,
		   double a[3],  double b[3], double c[3],
		   int mol_list[MAXATOM][MAXATOM],
		   int natom,
		   int nmolecule,
		   int *type,
		   char **element,
		   int bond_list[MAXATOM][MAXATOM],
                   int bond_time[MAXELE][MAXELE]) 
/** Recursively search for all atoms bonded to j, then all atoms 
    bonded to atoms bonded to j, etc. **/
{
  int k, l ;
  for ( k = 0 ; k < natom ; k++ ) {
    if ( is_bonded(j,k,x,y,z,a,b,c,type,element,bond_time) ) {
      bond_list[j][k] = 1 ;
      /**
      printf("Atoms %d and %d are bonded\n", j, k ) ;
      **/
      for ( l = 0 ; mol_list[nmolecule][l] >= 0 ; l++ ) {
	if ( mol_list[nmolecule][l] == k ) break ;
      }
      if ( mol_list[nmolecule][l] < 0 ) {
	/** We found a new entry **/
	mol_list[nmolecule][l] = k ;
	find_molecule(k,x,y,z,a,b,c,mol_list,natom,nmolecule,
		      type,element,bond_list,bond_time) ;
      }
    }
  }
}


int is_bonded(int j, int k, double *x, double *y, double *z,
	      double a[3], double b[3], double c[3],
	      int *type, char **element, int bond_time[MAXELE][MAXELE])
     /** Return 1 if j and k are bonded, 0 if not **/
     /** x, y, and z contain the atomic coordinates. **/
     /** a, b, and c are the rectangular box dimensions **/
{
  double r2 ;
  int i , btime , iele, jele ;

  if ( j == k ) return 0 ;

  r2 = wrap_atom(j,k, x, y, z, a, b, c, MAXWRAP, 0) ;
#if(0)
  printf("Testing Atoms %d (%s) and %d (%s) for bonding: r = %f\n", 
	 j, element[type[j]-1], k, element[type[k]-1], sqrt(r2) ) ;  
#endif
  /** The present passes the bonding criterion, now look at the past. */
  if ( r2 < r2bond(element[type[j]-1], element[type[k]-1]) ) {
    iele = ele_index(element[type[j]-1]);
    jele = ele_index(element[type[k]-1]);
    btime = bond_time[iele][jele] ;
    for ( i = 0 ; i < btime ; i++ ) {
      if ( xold[i][j] != NOTSET ) {
	/* Failing a bond test at any point in the past means that we are
	   not bonded. */
	r2 = wrap_atom(j,k, xold[i], yold[i], zold[i], a, b, c, MAXWRAP, 0) ;
	if ( r2 >= r2bond(element[type[j]-1], element[type[k]-1]) ) {
	  return(0) ;
	}
      }
    }
    return (1) ;
  }
  return(0) ;
}

  
  
double wrap_atom(int j, int k, double *x, double *y, double *z,
		 double a[3], double b[3], double c[3], int maxwrap,
		 int do_wrap)
     /** Wrap atom k so to be as close as possible to atom j if do_wrap is 
	 non-zero.
	 Returns the closest distance**2 found.  **/
{
  int l, m, n ;
  double x1[3] ;
  double bestdist = 0.0 , dist ;
  int bestl = 0 , bestm = 0 , bestn = 0 ;

  bestdist = 1.0e20 ;

  for ( l = -maxwrap ; l <= maxwrap ; l++ ) {
    for ( m = -maxwrap ; m <= maxwrap ; m++ ) {
      for ( n = -maxwrap ; n <= maxwrap ; n++ ) {

	x1[0] = x[j] + l * a[0] + m * b[0] + n * c[0] ;
	x1[1] = y[j] + l * a[1] + m * b[1] + n * c[1] ;
	x1[2] = z[j] + l * a[2] + m * b[2] + n * c[2] ;
	
	dist = (x1[0]-x[k]) * (x1[0] - x[k]) +
	  (x1[1] - y[k]) * (x1[1] - y[k]) +
	  (x1[2] - z[k]) * (x1[2] - z[k]) ;
	if ( dist < bestdist ) {
	  bestdist = dist ;
	  bestl = l ;
	  bestm = m ;
	  bestn = n ;
	}
      }
    }
  }
  if ( do_wrap ) {
#if(0)
    printf("bestdist = %f\n", bestdist) ;
    printf("bestl = %d, bestm = %d, bestn = %d\n", bestl, bestm, bestn) ;
    printf("dx = %f\n", bestl * a[0] + bestm * b[0] + bestn * c[0] ) ;
    printf("dy = %f\n", bestl * a[1] + bestm * b[1] + bestn * c[1] ) ;
    printf("dz = %f\n", bestl * a[2] + bestm * b[2] + bestn * c[2] ) ;
#endif
    x[k] -= bestl * a[0] + bestm * b[0] + bestn * c[0] ;
    y[k] -= bestl * a[1] + bestm * b[1] + bestn * c[1] ;
    z[k] -= bestl * a[2] + bestm * b[2] + bestn * c[2] ;
  }
  return bestdist ;
}



void wrap_com(double *x, double *y, double *z, 
	      double a[3], double b[3], double c[3],
	      int mol_list[MAXATOM][MAXATOM], int nmolecule, int natom,
	      double rcomx[MAXATOM], double rcomy[MAXATOM], 
	      double rcomz[MAXATOM], int type[MAXATOM],
	      char *element[MAXELE])
     /** Wrap each molecule so that it's center of mass (ignoring
	 atomic weight) is inside the box **/
{
  double rx, ry, rz, dx, dy, dz, rx0, ry0, rz0, mtot ;
  double mass ;
  double x1[2], y1[2], z1[2] ;
  double mass_cell ;

  int j , k, kk ;
  
  rx0 = 0.0 ;
  ry0 = 0.0 ;
  rz0 = 0.0 ;
  mass_cell = 0.0 ;
  for ( j = 0 ; j < nmolecule ; j++ ) {
    rx = 0.0 ;
    ry = 0.0 ;
    rz = 0.0 ;
    mtot = 0.0 ;
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      kk = mol_list[j][k] ;
      mass = atm_mass(element[type[mol_list[j][k]]-1]) ;
      mtot += mass ;
      rx += x[mol_list[j][k]] * mass ;
      ry += y[mol_list[j][k]] * mass ;
      rz += z[mol_list[j][k]] * mass ;
      /**
      printf("XX = %f YY = %f ZZ = %f Mass = %f type = %d\n", x[kk],
	     y[kk], z[kk], mass, type[mol_list[j][k]] ) ;
      printf("rx = %f rx = %f rz = %f\n", rx, ry, rz) ;
      **/
    }
    /* rx, ry, rz is the molecular center of mass. */
    rx /= mtot ;
    ry /= mtot ;
    rz /= mtot ;
    
    rcomx[j] = rx ;
    rcomy[j] = ry ;
    rcomz[j] = rz ;

    x1[0] = rx ;
    y1[0] = ry ;
    z1[0] = rz ;

    /** Calculate the center of the unit cell, with crystal coordinates
	(0.5, 0.5, 0.5) */
    x1[1] = 0.5 * a[0] + 0.5 * b[0] + 0.5 * c[0] ;
    y1[1] = 0.5 * a[1] + 0.5 * b[1] + 0.5 * c[1] ;
    z1[1] = 0.5 * a[2] + 0.5 * b[2] + 0.5 * c[2] ;
    
    /* Find the wrapped coordinates closest to the center of the unit cell. */
    wrap_atom(1,0,x1,y1,z1,a,b,c,MAXWRAP,1) ;

    dx = rx - x1[0] ;
    dy = ry - y1[0] ;
    dz = rz - z1[0] ;
#if(1)
    rcomx[j] -= dx ;
    rcomy[j] -= dy ;
    rcomz[j] -= dz ;

    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      x[mol_list[j][k]] -= dx ;
      y[mol_list[j][k]] -= dy ;
      z[mol_list[j][k]] -= dz ;
    }
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      mass = atm_mass(element[type[mol_list[j][k]]-1]) ;
      rx0 += x[mol_list[j][k]] * mass ;
      ry0 += y[mol_list[j][k]] * mass ;
      rz0 += z[mol_list[j][k]] * mass ;
      mass_cell += mass ;
    }
#endif
  }
  dx = rx0 / mass_cell ;
  dy = ry0 / mass_cell ;
  dz = rz0 / mass_cell ;
  for ( j = 0 ; j < nmolecule ; j++ ) {
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      x[mol_list[j][k]] -= dx ;
      y[mol_list[j][k]] -= dy ;
      z[mol_list[j][k]] -= dz ;
    }
  }
}

void wrap_molecules(double *x, double *y, double *z, 
		    double a[3], double b[3], double c[3],
		    int mol_list[MAXATOM][MAXATOM], int nmolecule,
		    int bond_list[MAXATOM][MAXATOM], int natom )
     /** Wrap the atoms of each molecule to be as close as possible
	 to bonded partners. **/
{
  double rx, ry, rz, dx, dy, dz, rx0, ry0, rz0, mtot ;
  double mass ;
  double x1[2], y1[2], z1[2] ;
  int wrapped[MAXATOM]  ;
  int j , k, kk, l, ll ;
  
  rx0 = 0.0 ;
  ry0 = 0.0 ;
  rz0 = 0.0 ;
  for ( j = 0 ; j < natom ; j++ ) {
    wrapped[j] = 0 ;
  }
  for ( j = 0 ; j < nmolecule ; j++ ) {
    rx = 0.0 ;
    ry = 0.0 ;
    rz = 0.0 ;
    mtot = 0.0 ;
    for ( k = 0 ; mol_list[j][k] >= 0 ; k++ ) {
      if ( wrapped[mol_list[j][k]] == 0 ) {
	wrap_pairs(x,y,z,a,b,c,mol_list,nmolecule,bond_list,natom,j,k, wrapped) ;
      }
    }
  }
}

void wrap_pairs(double *x, double *y, double *z, 
		double a[3], double b[3], double c[3],
		int mol_list[MAXATOM][MAXATOM], int nmolecule,
		int bond_list[MAXATOM][MAXATOM], int natom, int j, int k,
		int wrapped[MAXATOM])
     /** Wrap atoms bonded to the kth atom of the jth molecule. */
{  

  double rx, ry, rz, dx, dy, dz, rx0, ry0, rz0, mtot ;
  double mass ;
  double x1[2], y1[2], z1[2] ;
  int kk, l, ll ;

  kk = mol_list[j][k] ;
  for ( l = 0 ; mol_list[j][l] >= 0 ; l++ ) {
    ll = mol_list[j][l] ;
    if ( kk >= natom || ll >= natom ) {
      printf("Error in wrap_molecules: bad index\n") ;
      exit(1) ;
    }
    if ( kk != ll && bond_list[kk][ll] == 1 && wrapped[ll] == 0 ) {
      wrap_atom(kk,ll,x,y,z,a,b,c,MAXWRAP,1) ;
      wrapped[ll] = 1 ;
      wrap_pairs(x,y,z,a,b,c,mol_list,nmolecule,bond_list,natom,j,l,wrapped) ;
    }
  }
}

double r2bond(char *ea, char *eb)
     /** Return the square of the bond cutoff distance for 
	 element pair ea and eb.  This saves having to take the
     square root. **/
{
  int i, j ;

  i = ele_index(ea) ;
  j = ele_index(eb) ;
  if ( rbond[i][j] < 0.0 ) {
    printf("Error: bond %s-%s was not specified\n", ea, eb) ;
    exit(1) ;
  }

  return(rbond[i][j] * rbond[i][j]) ;
}

  
 
int nint(double num)
/* returns the nearest int to the double num */
{
 
    if (  fabs( num - (int ) num ) <  0.5)
        return( (int ) num ) ;
 
    else if ( num > 0 )
         return( (int ) num + 1 ) ;
    else
         return( (int ) num - 1 ) ;
}


void print_molecule(int mol_list[MAXATOM][MAXATOM],
		    int bond_list[MAXATOM][MAXATOM],
		    int *type, char **element, int mol, 
		    int read_charge,
		    double q[MAXATOM],
		    double x[MAXATOM], double y[MAXATOM], double z[MAXATOM],
		    double rcomx[MAXATOM], double rcomy[MAXATOM],
		    double rcomz[MAXATOM], int printed_atom[MAXATOM],
		    int dump_structures) 
{
  int k, j, jj, kk  ;
  int ele_count[MAXELE], bond_count[MAXELE][MAXELE] ;
  double charge = 0.0 ;
  double dipx, dipy, dipz ;
  double dipole, q0;
  char *ele ;
  char molecule_name[BUFSIZE] ;
  char buf[BUFSIZE] ;
  int natm_mol ;
  FILE *fmol ;

  memset(molecule_name, 0, BUFSIZE) ;
  memset(buf, 0, BUFSIZE) ;

  if ( dump_structures ) {
    strcat(molecule_name, "molecules/") ;
  }

  /** Calculate the overall stoichiometry **/
  for ( k = 0 ; k < MAXELE ; k++ ) {
    ele_count[k] = 0 ; 
  }
  dipx = 0.0 ;
  dipy = 0.0 ;
  dipz = 0.0 ;
  /**
  printf("COM: %f %f %f\n", rcomx[mol], rcomy[mol], rcomz[mol]) ;
  **/
  natm_mol = 0 ;
  for ( k = 0 ; mol_list[mol][k] >= 0 ; k++ ) {
    
    kk = mol_list[mol][k] ;
    /**
    printf("POS: %d %f %f %f\n", kk, x[kk], y[kk], z[kk]) ;
    **/
    q0 = - q[kk] ;
    ele = element[type[mol_list[mol][k]]-1] ;
    if ( strcmp(ele, "C") ==0 ) {
      q0 += 4 ;
    } else if ( strcmp(ele, "N") ==0 ) {
      q0 += 5 ;
    } else if ( strcmp(ele, "O") ==0) {
      q0 += 6 ;
    } else if ( strcmp(ele, "F") ==0) {
      q0 += 7 ;
    } else if ( strcmp(ele, "H") ==0) {
      q0 += 1 ;
    } else if ( strcmp(ele, "Cl") ==0) {
      q0 += 17 ;
    } else {
      printf("Error: charge for element %s is not known\n", ele) ;
    }
    charge += q0 ;
    {
      double ddx, ddy, ddz ;
      ddx = q0 * (x[kk]-rcomx[mol]);
      ddy = q0 * (y[kk]-rcomy[mol]) ;
      ddz = q0 * (z[kk]-rcomz[mol]) ;
      /* printf("z = %e, rcomz = %e q = %e\n", z[kk], rcomz[mol], q0) ; */
      dipx += ddx ;
      dipy += ddy ;
      dipz += ddz ;
      /*      printf("The dipole contribution for atom %d = %e %e %e\n", k, ddx, ddy, ddz) ; */
    }
    ele_count[type[mol_list[mol][k]]-1]++ ;
    natm_mol ++ ;
  }
  /* printf("The number of atoms in the molecule = %d\n", natm_mol) ; */

  for ( k = 0 ; mol_list[mol][k] >= 0 ; k++ ) {
    kk = mol_list[mol][k] ;
    ele = element[type[mol_list[mol][k]]-1] ;
    if ( strcmp(ele, "C") ==0 ) {
      q0 += 4 ;
    } else if ( strcmp(ele, "N") ==0 ) {
      q0 += 5 ;
    } else if ( strcmp(ele, "O") ==0 ) {
      q0 += 6 ;
    } else if ( strcmp(ele, "F") ==0 ) {
      q0 += 7 ;
    } else if ( strcmp(ele, "H") ==0 ) {
      q0 += 1 ;
    } else if ( strcmp(ele, "Cl") ==0 ) {
      q0 += 17 ;
    } else {
      printf("Error: charge for element %s is not known\n", ele) ;
    }
  }
  printf("Name: ") ;
  for ( k = 0 ; k < MAXELE ; k++ ) {
    if ( ele_count[k] > 0 ) {
      printf("%s%d ", element[k], ele_count[k]) ;
      sprintf(buf, "%s%d-", element[k], ele_count[k]) ;
      strncat(molecule_name, buf, BUFSIZE) ;
    }
  }

  /** Calculate the number of bonds **/
  for ( k = 0 ; k < MAXELE ; k++ ) {
    for ( j = 0 ; j < MAXELE ; j++ ) {
      bond_count[j][k] = 0 ;
    }
  }
  for ( k = 0 ; mol_list[mol][k] >= 0 ; k++ ) {
    kk = mol_list[mol][k] ;
    for ( j = 0 ; mol_list[mol][j] >= 0 ; j++ ) {
      jj = mol_list[mol][j] ;
      if ( bond_list[kk][jj] == 1 && jj > kk ) {
#if(0)	
	printf("Bonded: %d-%s %d-%s \n", kk, element[type[kk]-1],
	       jj, element[type[jj]-1]) ;
#endif
	if ( type[kk] > type[jj] ) {
	  bond_count[type[kk]-1][type[jj]-1]++ ;
	} else {
	  bond_count[type[jj]-1][type[kk]-1]++ ;
	}
      }
    }
  }

  for ( k = 0 ; k < MAXELE ; k++ ) {
    for ( j = 0 ; j < MAXELE ; j++ ) {
      if ( bond_count[j][k] > 0 ) {
	printf("%d(%s-%s) ", bond_count[j][k],
	       element[j], element[k]) ;
	sprintf(buf,"%d(%s-%s)-", bond_count[j][k],
	       element[j], element[k]) ;
	strncat(molecule_name, buf, BUFSIZE) ;
      }
    }
  }
  printf("\n") ;


  
  for ( k = 0 ; mol_list[mol][k] >= 0 ; ) {
    printf("   Atom list: ") ;
    for ( j = 0 ; j < 15 ; j++ ) {
      printf("%3d ", mol_list[mol][k] ) ;
      printed_atom[mol_list[mol][k]] = 1 ;
      k++ ;
      if ( mol_list[mol][k] < 0 ) {
	break ;
      }
    }
    printf("\n") ;
  }

  if ( dump_structures ) {
    strncat(molecule_name, ".xyz", BUFSIZE) ;
    fmol = fopen(molecule_name, "a") ;
    if ( fmol == NULL ) {
      printf("Error: could not open file %s\n", molecule_name) ;
      exit(1) ;
    }
    fprintf(fmol,"%d\n\n", k) ;
    for ( k = 0 ; mol_list[mol][k] >= 0 ; k++ ) {
      j = mol_list[mol][k] ;
      ele = element[type[mol_list[mol][k]]-1] ;
      fprintf(fmol,"%s %11.4f %11.4f %11.4f\n", ele, x[j], y[j], z[j]) ;
    }
    fclose(fmol) ;
  }

  /** 4.8 is a unit conversion factor. */
  dipole = 4.8 * sqrt(dipx * dipx + dipy * dipy + dipz * dipz) ;
  if ( read_charge ) {
    printf("   Charge = %f e\n", charge) ;
    printf("   Dipole moment = %f D\n", dipole) ;
  }

}
    

static void wrap_in_box(double a[3], double b[3], double c[3],
			double invbox[3][3],
			int natom, double x[MAXATOM], double y[MAXATOM], 
			double z[MAXATOM]) 
     /** This function wraps all the atoms into the primitive simulation
	 box. **/
{
  int j ;
  double ca, cb, cc ;

  for ( j = 0 ; j < natom ; j++ ) {
    ca = invbox[0][0] * x[j] + invbox[1][0] * y[j] + invbox[2][0] * z[j] ;
    cb = invbox[0][1] * x[j] + invbox[1][1] * y[j] + invbox[2][1] * z[j] ;
    cc = invbox[0][2] * x[j] + invbox[1][2] * y[j] + invbox[2][2] * z[j] ;

    ca -= nint(ca) ;
    cb -= nint(cb) ;
    cc -= nint(cc) ;

    x[j] = ca * a[0] + cb * b[0] + cc * c[0] ;
    y[j] = ca * a[1] + cb * b[1] + cc * c[1] ;
    z[j] = ca * a[2] + cb * b[2] + cc * c[2] ;

  }
}

void inversebox(double a[3], double b[3], double c[3],
  double invbox[3][3]) 
     /** Calculate the inverse box maxtrix, used for finding
	 crystal coordinates. */
{
  double xhlp ;
  xhlp=-a[2]*b[1]*c[0]+ a[1]*b[2]*c[0] ;
  xhlp=xhlp+a[2]*b[0]*c[1] ;
  xhlp=xhlp-a[0]*b[2]*c[1] ;
  xhlp=xhlp-a[1]*b[0]*c[2] ;
  xhlp=xhlp+a[0]*b[1]*c[2] ;

  invbox[0][0]=(-b[2]*c[1]+b[1]*c[2])/xhlp ;
  invbox[1][0]=(b[2]*c[0]-b[0]*c[2])/xhlp ;
  invbox[2][0]=(-b[1]*c[0]+b[0]*c[1])/xhlp ;

  invbox[0][1]=(a[2]*c[1]-a[1]*c[2])/xhlp ;
  invbox[1][1]=(-a[2]*c[0]+a[0]*c[2])/xhlp ;
  invbox[2][1]=(a[1]*c[0]-a[0]*c[1])/xhlp ;

  invbox[0][2]=(-a[2]*b[1]+a[1]*b[2])/xhlp ;
  invbox[1][2]=(a[2]*b[0]-a[0]*b[2])/xhlp ;
  invbox[2][2]=(-a[1]*b[0]+a[0]*b[1])/xhlp ;
}


double atm_mass(char *ea)
     /* Return the atomic mass of the given element. */
{
  if ( strcmp(ea, "H") ==0 ) {
    return(1.0) ;
  } else if ( strcmp(ea, "C") ==0 ) {
    return(12.0) ;
  } else if ( strcmp(ea, "N") ==0 ) {
    return(14.0) ;
  } else if ( strcmp(ea, "O") ==0 ) {
    return(16.0) ;
  } else if ( strcmp(ea, "F") ==0 ) {
    return(19.0) ;
  } else if ( strcmp(ea, "Cl") ==0 ) {
    return(35.5) ;
  } else {
    printf("Error: element %s is unknown\n", ea ) ;
    exit(1) ;
  }
  return(1.0) ;
}

void read_bonds(double rbond[MAXELE][MAXELE], int *duration, int *xyz_copies, double *time_step,
		int *dump_structures, int bond_time[MAXELE][MAXELE])
     /** Read in bond distances and other parameters.
	 rbond - The maximum bond distance.
	 duration - The bond lifetime in simulation steps.
	 xyz_copies - The number of copies in each direction for the xyz output file.
	 time_step - The time between each frame stored.
	 dump_structures - Whether to write out individual molecule structure files.
     **/
{
  FILE *fbond, *ftime ;
  char *ele1, *ele2 ;
  char buf1[BUFSIZ], buf2[BUFSIZ] ;
  double r ;
  int i, j, k, btime ;

  for ( i = 0 ; i < MAXELE ; i++ ) {
    for ( j = 0 ; j < MAXELE ; j++ ) {
      /** Put in a negative value so that unspecified bonds can be caught. */
      rbond[i][j] = -1.0 ;
      bond_time[i][j] = -1 ;
    }
  }

  fbond = fopen("bonds.dat","r") ;
  if ( fbond == NULL ) {
    printf("Error: could not open bonds.dat\n") ;
    exit(1) ;
  }
  fscanf(fbond, "%d", duration) ;
  if (*duration >= 0) {
    printf("Bonds require a lifetime of %d steps\n", *duration)  ;
  } else {
    ftime = fopen("bond_times.dat","r") ;
    if (ftime == NULL) {
      printf ("Error: could not open bond_times.dat\n");
      exit(1);
    }
  }
  fscanf(fbond, "%d", xyz_copies) ;
  printf("%d replicas of the simulation cell will be propagated in each direction\n",
	 *xyz_copies) ;

  fscanf(fbond, "%lf\n", time_step) ;
  printf("The time between frames read in = %11.7e seconds\n", *time_step) ;
  if ( *time_step < 0.0 || *time_step > 1.0e-10 ) {
    printf("The time step was given as %21.14e seconds.  This seems strange\n") ;
    exit(1) ;
  }

  fscanf(fbond, "%d", dump_structures) ;
  if ( *dump_structures != 0 ) {
    printf("Individual molecule xyz files will be created in the molecules directory.\n") ;
  } else {
    printf("Individual molecule xyz file will not be created.\n") ;
  }
  printf("Bond distances:\n") ;

  do {
    if ( fscanf(fbond, "%s %s %lf", buf1, buf2, &r) != 3 ) {
      break ;
    }
    ele1 = buf1 ;
    ele2 = buf2 ;
    i = ele_index(ele1) ;
    j = ele_index(ele2) ;
    if ( j < i ) {
      k = i ;
      i = j ;
      j = k ;
    }
    printf("%s %s %f\n", ele1, ele2, r) ;
    if ( rbond[i][j] > 0.0 || rbond[j][i] > 0.0 ) {
      printf("Error: Bond distance for %s %s has already been specified\n",
	     ele1, ele2) ;
      exit(1) ;
    }
    rbond[i][j] = r ;
    if (*duration >=0 ) {
      bond_time[i][j] = *duration ;
    }
    
  } while ( ! feof(fbond) && ! ferror(fbond)) ;
  fclose(fbond) ;


  printf("Bond time cutoffs:\n") ;
  if (*duration >= 0) {
    for (i = 0; i < MAXELE; i++) {
      for (j = 0; j < MAXELE; j++) {
        if (i == 0) {
          ele1 = "H";
        } else if ( i == 1) {
          ele1 = "C";
        } else if ( i == 2) {
          ele1 = "N";
        } else if ( i == 3) {
          ele1 = "O";
        } else if ( i == 4) {
          ele1 = "F";
        } else if ( i == 5) {
          ele1 = "Cl";
        }
        if (j == 0) {
          ele2 = "H";
        } else if ( j == 1) {
          ele2 = "C";
        } else if ( j == 2) {
          ele2 = "N";
        } else if ( j == 3) {
          ele2 = "O";
        } else if ( j == 4) {
          ele2 = "F";
        } else if ( j == 5) {
          ele2 = "Cl";
        }
        if (bond_time[i][j] != -1) {
          printf("%s %s %d\n", ele1, ele2, bond_time[i][j]) ;
        }
      }
    }
  } else {
    do {
      if ( fscanf(ftime, "%s %s %d", buf1, buf2, &btime) != 3 ) {
        break ;
      }
      ele1 = buf1 ;
      ele2 = buf2 ;
      i = ele_index(ele1) ;
      j = ele_index(ele2) ;
      if ( j < i ) {
        k = i ;
        i = j ;
        j = k ;
      }
      printf("%s %s %d\n", ele1, ele2, btime) ;
      if ( bond_time[i][j] > 0.0 || bond_time[j][i] > 0.0 ) {
        printf("Error: Bond time cutoff for %s %s has already been specified\n",
	       ele1, ele2) ;
        exit(1) ;
      }
      bond_time[i][j] = btime ;
    
    } while ( ! feof(ftime) && ! ferror(ftime)) ;
    fclose(ftime) ;
  }

  for ( i = 0 ; i < MAXELE ; i++ ) {
    for ( j = i + 1 ; j < MAXELE ; j++ ) {
      rbond[j][i] = rbond[i][j] ;
      bond_time[j][i] = bond_time[i][j] ;
    }
  }
}

int ele_index(char *ea)
     /** Returns the index of the given element */
{
  int i ;
  if ( strcmp(ea, "H") == 0 ) {
    i = 0 ;
  } else if ( strcmp(ea, "C") == 0 ) {
    i = 1 ;
  } else if (strcmp(ea, "N") == 0 ) {
    i = 2 ;
  } else if (strcmp(ea, "O") == 0 ) {
    i = 3 ; 
  } else if (strcmp(ea, "F") == 0 ) {
    i = 4 ; 
  } else if (strcmp(ea, "Cl") == 0 ) {
    i = 5 ; 
  } else {
    printf("Error: element %s is unknown\n", ea ) ;
    exit(1) ;
  }
  return i ;
}

double boxvol(double a[3], double b[3], double c[3])
     /* Calculates the volume of the simulation box, given the
	box cell vectors a, b, and c. */
{

  double bv ;
  double xhlp ;
  double ax, ay, az, bx, by, bz, cx, cy, cz ;

  ax = a[0] ;
  ay = a[1] ;
  az = a[2] ;
  bx = b[0] ;
  by = b[1] ;
  bz = b[2] ;
  cx = c[0] ;
  cy = c[1] ;
  cz = c[2] ;

  xhlp=-az*by*cx+ ay*bz*cx ;
  xhlp=xhlp+az*bx*cy ;
  xhlp=xhlp-ax*bz*cy ;
  xhlp=xhlp-ay*bx*cz ;
  xhlp=xhlp+ax*by*cz ;

  bv = fabs(xhlp) ;
  return bv;
}

