#!/bin/tcsh

foreach kpoint ( 04*k? )

cd $kpoint

  echo $kpoint

    awk '{print $2, $10}' total >> ../eig.dat

  cd ..

end
