#!/bin/tcsh

foreach machine (cl001 cl002 cl003 cl004 cl005 cl006)

  echo $machine
  ssh $machine /usr/local/bin/qclean-mx.sh

end
