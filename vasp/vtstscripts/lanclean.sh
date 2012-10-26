#!/bin/sh

  echo "   ****   "
  echo " Cleaning up after a lanczos vasp run "
  echo "   ****   " 

  Bin=`dirname "$0"`
  ZIP=${VTST_ZIP-gzip}

  rm -f DOSCAR EIGENVAL IBZKPT OSZICAR PCDAT vasprun.xml CHG

  if [ ${VTST_STDOUT} ] 
  then
    if [ -s ${VTST_STDOUT} ]
    then
      mv ${VTST_STDOUT} $1
    fi
  fi
  if [ ${VTST_STDERR} ]
  then
    rm -f ${VTST_STDERR}
  fi

  if [ -s WAVECAR ] 
  then
    mv WAVECAR $1 ; "$ZIP" $1/WAVECAR &
  else
    rm -f WAVECAR
  fi
  if [ -s CHGCAR ] 
  then
    mv CHGCAR $1 ; "$ZIP" $1/CHGCAR &
  else
    rm -f CHGCAR
  fi

  if [ -s DISPLACECAR ]
  then
    rm -f DISPLACECAR
  fi 

  if [ -s NEWMODECAR ]  
  then
    cp MODECAR NEWMODECAR $1
    mv NEWMODECAR MODECAR
  fi

  if [ -s lanczos.out ]      
  then
    mv lanczos.out $1 ; "$ZIP" $1/lanczos.out &
  else
    rm -f lanczos.out
  fi
 
  rm PI*
  cp POSCAR CONTCAR INCAR KPOINTS $1

  "$Bin/xdat2xyz.pl" > /dev/null ;
  if [ -s movie.xyz ]
  then
     mv movie.xyz $1
  fi
  "$Bin/vef.pl" > /dev/null ;
  if [ -s fe.dat ]
  then
     mv fe.dat $1
  fi
  if [ -s vaspout.eps ]
  then
     mv vaspout.eps $1
  fi
  if [ -s XDATCAR ]
  then
     mv XDATCAR $1 ; "$ZIP" $1/XDATCAR &
  fi
  if [ -s OUTCAR ]
  then
     mv OUTCAR $1 ; "$ZIP" $1/OUTCAR &
  fi

  "$Bin/vasp2con.pl" POSCAR > /dev/null
  "$Bin/con2xyz.pl" POSCAR.con > /dev/null 
  mv POSCAR.con POSCAR.xyz $1 ;
  "$Bin/vasp2con.pl" CONTCAR > /dev/null 
  "$Bin/con2xyz.pl" CONTCAR.con > /dev/null 
  mv CONTCAR.con CONTCAR.xyz $1 ;

  mv CONTCAR POSCAR
