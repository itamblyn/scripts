#!/bin/tcsh

set file=$1

setenv g09root /usr/local/gaussian
setenv GAUSS_SCRDIR /globalscratch/$USER/
source /usr/local/gaussian/g09/bsd/g09.login
formchk $file
