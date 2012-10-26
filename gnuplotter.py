#!/usr/bin/env python
"""
gnuplotter.sh
Brian Boates

To enable the use of wildcards (i.e. *.dat) in gnuplotting
"""
import os, sys, glob, commands


def main():

    # Get filenames from user
    try:
        fnames = sys.argv[1]
    except:
        print '\n usage: '+sys.argv[0]+' "*.dat"\n'
        sys.exit(1)

    # Write a generic gnuplot .plt file
    out = open('new.plt','w')

    out.write('#!/usr/bin/gnuplot -persist'+'\n')
    out.write('#'+'\n')
    out.write('#    '+'\n')
    out.write('#    	G N U P L O T'+'\n')
    out.write('#    	Version 4.0 patchlevel 0'+'\n')
    out.write('#    	last modified Thu Apr 15 14:44:22 CEST 2004'+'\n')
    out.write('#    	System: Linux 2.6.18-194.el5'+'\n')
    out.write('#    '+'\n')
    out.write('#    	Copyright (C) 1986 - 1993, 1998, 2004'+'\n')
    out.write('#    	Thomas Williams, Colin Kelley and many others'+'\n')
    out.write('#    '+'\n')
    out.write('#    	This is gnuplot version 4.0.  Please refer to the documentation'+'\n')
    out.write('#    	for command syntax changes.  The old syntax will be accepted'+'\n')
    out.write('#    	throughout the 4.0 series, but all save files use the new syntax.'+'\n')
    out.write('#    '+'\n')
    out.write('#    	Type `help` to access the on-line reference manual.'+'\n')
    out.write('#    	The gnuplot FAQ is available from'+'\n')
    out.write('#    		http://www.gnuplot.info/faq/'+'\n')
    out.write('#    '+'\n')
    out.write('#    	Send comments and requests for help to'+'\n')
    out.write('#    		<gnuplot-info@lists.sourceforge.net>'+'\n')
    out.write('#    	Send bugs, suggestions and mods to'+'\n')
    out.write('#    		<gnuplot-bugs@lists.sourceforge.net>'+'\n')
    out.write('#    '+'\n')
    out.write('# set terminal x11 '+'\n')
    out.write('# set output'+'\n')
    out.write('unset clip points'+'\n')
    out.write('set clip one'+'\n')
    out.write('unset clip two'+'\n')
    out.write('set bar 1.000000'+'\n')
    out.write('set border 31 lt -1 lw 1.000'+'\n')
    out.write('set xdata'+'\n')
    out.write('set ydata'+'\n')
    out.write('set zdata'+'\n')
    out.write('set x2data'+'\n')
    out.write('set y2data'+'\n')
    out.write('set timefmt x "%d/%m/%y,%H:%M"'+'\n')
    out.write('set timefmt y "%d/%m/%y,%H:%M"'+'\n')
    out.write('set timefmt z "%d/%m/%y,%H:%M"'+'\n')
    out.write('set timefmt x2 "%d/%m/%y,%H:%M"'+'\n')
    out.write('set timefmt y2 "%d/%m/%y,%H:%M"'+'\n')
    out.write('set timefmt cb "%d/%m/%y,%H:%M"'+'\n')
    out.write('set boxwidth'+'\n')
    out.write('set style fill empty border'+'\n')
    out.write('set dummy x,y'+'\n')
    out.write('set format x "% g"'+'\n')
    out.write('set format y "% g"'+'\n')
    out.write('set format x2 "% g"'+'\n')
    out.write('set format y2 "% g"'+'\n')
    out.write('set format z "% g"'+'\n')
    out.write('set format cb "% g"'+'\n')
    out.write('set angles radians'+'\n')
    out.write('unset grid'+'\n')
    out.write('set key title ""'+'\n')
    out.write('set key right top Right noreverse enhanced box linetype -2 linewidth 1.000 samplen 4 spacing 1 width 0 height 0 autotitles'+'\n')
    out.write('unset label'+'\n')
    out.write('unset arrow'+'\n')
    out.write('unset style line'+'\n')
    out.write('unset style arrow'+'\n')
    out.write('unset logscale'+'\n')
    out.write('set offsets 0, 0, 0, 0'+'\n')
    out.write('set pointsize 1'+'\n')
    out.write('set encoding default'+'\n')
    out.write('unset polar'+'\n')
    out.write('unset parametric'+'\n')
    out.write('unset decimalsign'+'\n')
    out.write('set view 60, 30, 1, 1'+'\n')
    out.write('set samples 100, 100'+'\n')
    out.write('set isosamples 10, 10'+'\n')
    out.write('set surface'+'\n')
    out.write('unset contour'+'\n')
    out.write("set clabel '%8.3g'"+"\n")
    out.write('set mapping cartesian'+'\n')
    out.write('set datafile separator whitespace'+'\n')
    out.write('unset hidden3d'+'\n')
    out.write('set cntrparam order 4'+'\n')
    out.write('set cntrparam linear'+'\n')
    out.write('set cntrparam levels auto 5'+'\n')
    out.write('set cntrparam points 5'+'\n')
    out.write('set size ratio 0 1,1'+'\n')
    out.write('set origin 0,0'+'\n')
    out.write('set style data linespoints'+'\n')
    out.write('set style function lines'+'\n')
    out.write('set xzeroaxis lt -2 lw 1.000'+'\n')
    out.write('set yzeroaxis lt -2 lw 1.000'+'\n')
    out.write('set x2zeroaxis lt -2 lw 1.000'+'\n')
    out.write('set y2zeroaxis lt -2 lw 1.000'+'\n')
    out.write('set tics in'+'\n')
    out.write('set ticslevel 0.5'+'\n')
    out.write('set ticscale 1 0.5'+'\n')
    out.write('set mxtics default'+'\n')
    out.write('set mytics default'+'\n')
    out.write('set mztics default'+'\n')
    out.write('set mx2tics default'+'\n')
    out.write('set my2tics default'+'\n')
    out.write('set mcbtics default'+'\n')
    out.write('set xtics border mirror norotate autofreq '+'\n')
    out.write('set ytics border mirror norotate autofreq '+'\n')
    out.write('set ztics border nomirror norotate autofreq '+'\n')
    out.write('set nox2tics'+'\n')
    out.write('set noy2tics'+'\n')
    out.write('set cbtics border mirror norotate autofreq '+'\n')
    out.write('set title "" 0.000000,0.000000  font ""'+'\n')
    out.write('set timestamp "" bottom norotate 0.000000,0.000000  ""'+'\n')
    out.write('set rrange [ * : * ] noreverse nowriteback  # (currently [0.00000:10.0000] )'+'\n')
    out.write('set trange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )'+'\n')
    out.write('set urange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )'+'\n')
    out.write('set vrange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )'+'\n')
    out.write('set xlabel "" 0.000000,0.000000  font ""'+'\n')
    out.write('set x2label "" 0.000000,0.000000  font ""'+'\n')
    out.write('set xrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )'+'\n')
    out.write('set x2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )'+'\n')
    out.write('set ylabel "" 0.000000,0.000000  font ""'+'\n')
    out.write('set y2label "" 0.000000,0.000000  font ""'+'\n')
    out.write('set yrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )'+'\n')
    out.write('set y2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )'+'\n')
    out.write('set zlabel "" 0.000000,0.000000  font ""'+'\n')
    out.write('set zrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )'+'\n')
    out.write('set cblabel "" 0.000000,0.000000  font ""'+'\n')
    out.write('set cbrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )'+'\n')
    out.write('set zero 1e-08'+'\n')
    out.write('set lmargin -1'+'\n')
    out.write('set bmargin -1'+'\n')
    out.write('set rmargin -1'+'\n')
    out.write('set tmargin -1'+'\n')
    out.write('set locale "C"'+'\n')
    out.write('set pm3d scansautomatic flush begin noftriangles nohidden3d implicit corners2color mean'+'\n')
    out.write('unset pm3d'+'\n')
    out.write('set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB '+'\n')
    out.write('set palette rgbformulae 7, 5, 15'+'\n')
    out.write('set colorbox default'+'\n')
    out.write('set colorbox vertical origin 0.9,0.2 size 0.1,0.63 bdefault'+'\n')
    out.write('set loadpath '+'\n')
    out.write('set fontpath '+'\n')
    out.write('set fit noerrorvariables'+'\n')


    # Get the filenames the user asked for
    fs = glob.glob(fnames)
    fs.sort()

    # Create the plot line for gnuplot
    line = 'plot '
    for f in fs:
        line += '"'+f+'",'
    line = line[:-1] # remove last comma

    # write the plot line to gnuplot file
    out.write(line+'\n')

    out.write('#    EOF'+'\n')
    out.close()  # end of .plt file

    # print user instructions
    print 
    print "View your data by either entering gnuplot and typing 'load new.plt'"
    print
    print " === OR === by copying the below line directly into gnuplot yourself"
    print
    print line
    print


if __name__ == '__main__':
    main()
