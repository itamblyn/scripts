#!/usr/bin/python

import numpy

def read_eigenvalues(subtractFermi):

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

def main():

    subtractFermi = True

    eigenvalue_list, wkpt_array, nkpt, neigen_per_kpt = read_eigenvalues(subtractFermi) 

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
