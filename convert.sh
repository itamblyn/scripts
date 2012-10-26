#!/bin/bash

echo "18839" > coor.xyz
echo "110., 110., 1000.0" >> coor.xyz

awk '/O/{print "O ", $2, $3, $4}' $1 >> coor.xyz
awk '/SI/{print "Si ", $2, $3, $4}' $1 >> coor.xyz
awk '/H/{print "H ", $2, $3, $4}' $1 >> coor.xyz

echo "coor.xyz" > RDF.in
echo "RDF.O-Si.dat" >> RDF.in
echo "O" >> RDF.in
echo "Si" >> RDF.in
echo "10.0" >> RDF.in
echo "0.025" >> RDF.in
echo "1" >> RDF.in

./RDF < RDF.in 

echo "coor.xyz" > RDF.in
echo "RDF.O-O.dat" >> RDF.in
echo "O" >> RDF.in
echo "O" >> RDF.in
echo "10.0" >> RDF.in
echo "0.025" >> RDF.in
echo "1" >> RDF.in

./RDF < RDF.in 


echo "coor.xyz" > RDF.in
echo "RDF.Si-Si.dat" >> RDF.in
echo "Si" >> RDF.in
echo "Si" >> RDF.in
echo "10.0" >> RDF.in
echo "0.025" >> RDF.in
echo "1" >> RDF.in

./RDF < RDF.in 

echo "coor.xyz" > RDF.in
echo "RDF.O-H.dat" >> RDF.in
echo "O" >> RDF.in
echo "H" >> RDF.in
echo "10.0" >> RDF.in
echo "0.025" >> RDF.in
echo "1" >> RDF.in

./RDF < RDF.in 

