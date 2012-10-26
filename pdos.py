#!/usr/bin/env python

import sys
sys.path.append("/usr/local/numpy.1.0.4/lib64/python2.3/site-packages/")
import numpy
import commands

try:
    input_filename = sys.argv[1]
except IndexError:
    print '\nusage: ' + sys.argv[0] + ' input_filename'
    print '\nexiting...\n'
    sys.exit(0)

core_valence_divide = 401 # <- this corresponds to -20 eV 

command_line_counter = commands.getoutput('wc -l ' + input_filename).split()

if len(command_line_counter) != 2:
    print 'Error determining file size'
else:
    number_of_lines = int(command_line_counter[0])

outputFile_dos = open('dos.dat', 'w')
outputFile_pdos = open('pdos.dat', 'w')

inputFile = open(input_filename, 'r')

for i in range(5):
    inputFile.readline()

inline_total = inputFile.readline()

emin_total = float(inline_total.split()[1])
emax_total = float(inline_total.split()[0])
enum_total = int(inline_total.split()[2])
efermi_total = float(inline_total.split()[3])

dos = numpy.zeros((enum_total, 3), dtype = numpy.float) # E, dos, idos

for i in range(enum_total):
    inline_total = inputFile.readline().split()
    dos[i][0] = float(inline_total[0])  #- efermi_total 
#    dos[i][0] = float(inline_total[0]) - efermi_total 
    dos[i][1], dos[i][2] = float(inline_total[1]), float(inline_total[2])

for row in dos:
    for element in row:
        outputFile_dos.write(str(element) + ' ')
    outputFile_dos.write('\n')

############

inline_project = inputFile.readline()

emin_project = float(inline_project.split()[1])
enum_project = int(inline_project.split()[2])
efermi_project = float(inline_project.split()[3])

spacing = (emax_total - emin_total)/enum_total
efermi_bin = int((efermi_total - emin_total)/spacing)

pdos = numpy.zeros((enum_project, 1 + 3*2), dtype = numpy.float) # E, s, p, d

for i in range(len(pdos)):

   pdos[i][0] = dos[i][0]  # fills in the E column of pdos

for i in range(enum_project):
    inline_project = inputFile.readline().split()
    pdos[i][1] += float(inline_project[1])
    pdos[i][2] += float(inline_project[2])
    pdos[i][3] += float(inline_project[3])

for i in numpy.arange(1,enum_project):
    pdos[i][4] = pdos[i - 1][4] + pdos[i][1]
    pdos[i][5] = pdos[i - 1][5] + pdos[i][2]
    pdos[i][6] = pdos[i - 1][6] + pdos[i][3]


energies  = numpy.transpose(pdos)[0]
s_project = numpy.transpose(pdos)[1]
p_project = numpy.transpose(pdos)[2]
d_project = numpy.transpose(pdos)[3]

energies_core  =  energies[0:core_valence_divide]
s_project_core = s_project[0:core_valence_divide]
p_project_core = p_project[0:core_valence_divide]
d_project_core = d_project[0:core_valence_divide]


energies_valence  =  energies[core_valence_divide:efermi_bin]
s_project_valence = s_project[core_valence_divide:efermi_bin]
p_project_valence = p_project[core_valence_divide:efermi_bin]
d_project_valence = d_project[core_valence_divide:efermi_bin]

spd_sum_core    = sum(s_project_core)    + sum(p_project_core)    + sum(d_project_core)
spd_sum_valence = sum(s_project_valence) + sum(p_project_valence) + sum(d_project_valence)

s_fraction_core = sum(s_project_core)/spd_sum_core
p_fraction_core = sum(p_project_core)/spd_sum_core

s_fraction_valence = sum(s_project_valence)/spd_sum_valence
p_fraction_valence = sum(p_project_valence)/spd_sum_valence

print
print "%s core", s_fraction_core*100
print "%p core", p_fraction_core*100
print
print "%s valence", s_fraction_valence*100
print "%p valence", p_fraction_valence*100
print
print "efermi_total = ", efermi_total
print "spacing = ", spacing
print

for row in pdos:
    for element in row:
        outputFile_pdos.write(str(element) + ' ')
    outputFile_pdos.write('\n')
