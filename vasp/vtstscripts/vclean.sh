#!/bin/sh

  echo "   ****   "
  echo " Cleaning up after a normal vasp run "
  echo "   ****   " 

  Bin=`dirname "$0"`
  ZIP=${VTST_ZIP-gzip}

  rm -f EIGENVAL IBZKPT OSZICAR PCDAT vasprun.xml CHG

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
  if [ -s PROCAR ]
  then
    mv PROCAR $1 ; "$ZIP" $1/PROCAR &
  else
    rm -f PROCAR
  fi
  if [ -s DOSCAR ] 
  then
    mv DOSCAR $1 ; "$ZIP" $1/DOSCAR &
  else
    rm -f DOSCAR
  fi

  "$Bin/xdat2xyz.pl" > /dev/null ;
  if [ -s movie.xyz ]
  then
     mv movie.xyz $1
  fi
  "$Bin/vef.pl" > /dev/null ;
  if [ -s vaspout.eps ]
  then
     mv fe.dat vaspout.eps $1
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
  cp POSCAR CONTCAR INCAR KPOINTS $1

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

  mv CONTCAR POSCAR
