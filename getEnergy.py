#!/usr/bin/python

import sys, commands

combined_energy = float(commands.getoutput('grep "Convergence criterion met" ../output.out | tail -1 ').split()[1])
iso_counter = int(commands.getoutput('grep "Convergence criterion met" output.out | wc -l '))

iso_energy_line = commands.getoutput('grep "Convergence criterion met" output.out').split()

total = combined_energy



iso_energy = []

for i in range(iso_counter):

    iso_energy.append(float(iso_energy_line[1 + i*6]))
    total -= iso_energy[-1]

print combined_energy, iso_energy
print -total*627.509469


