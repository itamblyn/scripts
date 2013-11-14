#!/bin/bash

etot=`grep Conv output.out | head -1 | awk '{print $2}'`

awk -v e=$etot '/Conv/{print ($2 - e)*27.2}' output.out
