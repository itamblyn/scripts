#!/bin/sh
#
# preparation of the mpi machine file
#
# usage: startmpi.sh <pe_hostfile>
#

PeHostfile2MachineFile()
{
   cat $1 | while read line; do
      host=`echo $line|cut -f1 -d" "|cut -f1 -d"."`
      nslots=`echo $line|cut -f2 -d" "`
      i=1
      while [ $i -le $nslots ]; do
         echo $host
         i=`expr $i + 1`
      done
   done
}

#
# on success the job will find a machine-file in $TMPDIR/machines
# 

me=`basename $0`

# test number of args
if [ $# -ne 1 ]; then
   echo "$me: got wrong number of arguments" >&2
   exit 1
fi

# get arguments
pe_hostfile=$1

# ensure pe_hostfile is readable
if [ ! -r $pe_hostfile ]; then
   echo "$me: can't read $pe_hostfile" >&2
   exit 1
fi

machines="$TMPDIR/machines"
PeHostfile2MachineFile $pe_hostfile >> $machines

# signal success to caller
exit 0
