#!/usr/bin/python

import numpy, commands

def read_vasp_eigenvalues(subtractFermi):

    inputFile_EIGENVAL = open('EIGENVAL', 'r')
    inputFile_IBZKPT = open('IBZKPT', 'r')
    inputFile_DOSCAR = open('DOSCAR', 'r')

    for i in range(5):
        inputFile_EIGENVAL.readline()
        inputFile_DOSCAR.readline()
    for i in range(3):
        inputFile_IBZKPT.readline()

    efermi = float(inputFile_DOSCAR.readline().split()[3])

    line = inputFile_EIGENVAL.readline()

    nelectrons     = int(line.split()[0])
    nkpt           = int(line.split()[1])
    neigen_per_kpt = int(line.split()[2])

    print nelectrons, ' electrons'
    print nkpt, ' kpoints'
    print neigen_per_kpt, ' eigenvalues per kpoint'

    print 'Fermi level at: ', efermi

    wkpt_array = numpy.zeros(nkpt, dtype=float)
    eigenvalue_array = []

    for i in range(nkpt):

        eigenvalue_array.append([])
        inputFile_EIGENVAL.readline()   # skips line before data
        wkpt = float(inputFile_EIGENVAL.readline().split()[3])
        wkpt_array[i] = wkpt

        for j in range(neigen_per_kpt):
            eigenvalue = float(inputFile_EIGENVAL.readline().split()[1])
            if subtractFermi == True: eigenvalue_array[-1].append(eigenvalue - efermi)
            else: eigenvalue_array[-1].append(eigenvalue)

    eigenvalue_list = []

    for i in range(nkpt):
        for eigenvalue in eigenvalue_array[i]:
            eigenvalue_list.append(eigenvalue)

    return eigenvalue_list, wkpt_array, nkpt, neigen_per_kpt

def read_paratec_eigenvalues(subtractFermi):

    efermi = float(commands.getoutput('awk \'/Fermi/{print $7}\' OUT | tail -1'))
    print 'efermi = ', efermi 
    inputFile = open('SCF_KPOINTS')

    nkpt = int(inputFile.readline().split()[0])  
    
    wkpt_array = numpy.zeros(nkpt, dtype=float)
    eigenvalue_list = []

    neigen_per_kpt= int(commands.getoutput('grep "m\=" OUT | head -1 ').split()[3])
    print 'neigen_per_kpt =', neigen_per_kpt

    for i in range(nkpt):
        wkpt_array[i] = float(inputFile.readline().split()[4])
 
    inputFile.close()
    grep_lines = neigen_per_kpt*2/7 + 1 # two lines per, 7 columns, and one line which says kpoint
    
    command_line_counter = commands.getoutput('grep -A ' + str(grep_lines) + ' k-point OUT | tac | grep "\." | grep -v "(" > grepfile')

    inputFile = open('grepfile')

    # read in ALL eigenvalues (it will be easier this way)

    tmp_array = []

    for line in inputFile.readlines():
        for element in line.split():
            tmp_array.append(float(element))

    counter = 0
    for inkpt in range(nkpt):
        eig_array = []
        for ineigen_per_kpt in range(neigen_per_kpt):
            eig_array.append(tmp_array[counter])
            counter +=1

        eig_array.sort()
        eigenvalue_list = eig_array + eigenvalue_list


    eigenvalue_list = numpy.array(eigenvalue_list)
    if subtractFermi == True: eigenvalue_list -= efermi

    return eigenvalue_list, wkpt_array, nkpt, neigen_per_kpt

def main():

    subtractFermi = True
    detect_paratec = int(commands.getoutput('grep paratecSGL OUT | wc -l '))
    
    if detect_paratec == 1:
        print 'Detected paratec'
        eigenvalue_list, wkpt_array, nkpt, neigen_per_kpt = read_paratec_eigenvalues(subtractFermi) 
    else:
        print 'Detected vasp'
        eigenvalue_list, wkpt_array, nkpt, neigen_per_kpt = read_vasp_eigenvalues(subtractFermi) 
 
    command = commands.getoutput('rm -f grepfile')
  
    print len(eigenvalue_list), ' eigenvalues were read'

#    g = Gnuplot.Gnuplot()
#    g('set data style linespoints') 
#    g.plot(sorted_array)

#    raw_input('Please press return to continue...\n')

    outputFile = open('eig.dat', 'w')

    outputFile.write('# nkptgw: ' + str(nkpt) + ' neig: ' + str(neigen_per_kpt) + '\n')
    outputFile.write('# E(DFT) pulse\n')

    for element in eigenvalue_list:
        outputFile.write(str(element) + ' ' + str(1.0) + ' \n')

    outputFile.close()

    outputFile = open('wtk.dat', 'w')

    i = 0
    for element in wkpt_array:
        i += 1
        outputFile.write(str(element) + ' ')
        if ( i % 6 == 0 ): outputFile.write('\n')

    outputFile.write('\n')

    outputFile.close()

    if subtractFermi == True: print 'Note, the fermi level has been subtracted'

if __name__ == '__main__':
   main()
