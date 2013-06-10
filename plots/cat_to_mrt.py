#! /usr/bin/env python
"""
python script to parse output of avg_spectrum.py into input suitable for the MRT creator at
http://authortools.aas.org/MRT/upload.html

input a catalog using usual aipy convention for positions of the sources
"""
import numpy as n,sys,optparse,aipy as a,ephem
from pylab import *


o = optparse.OptionParser()
o.set_usage('cat_to_mrt.py [options] psa64_pic_stripe.txt')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, src=True)
o.add_option('--filename',default='raw_table.txt',
    help='output MRT file. Default=raw_table.txt')
opts,args = o.parse_args(sys.argv[1:])

lines = open(args[-1]).readlines()
print 
freqs = map(float,lines[0].split('=')[1].split(','))

lines = [line.split() for line in lines if not line.startswith('#')]
srclist = [l[0] for l in lines]

#if not opts.cat is None:
srclist,cutoff,catalogs = a.scripting.parse_srcs(','.join(srclist), opts.cat)
if not opts.cal is None:
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
else:
    cat = a.src.get_catalog(srclist, cutoff, catalogs)
o = a.scripting.get_null_aa()
cat.compute(o)

#first output the MRT input table
#MRTfilename = args[-1][:-4]+'_MRT.txt'
            
F = open(opts.filename,'w')
recs = []
picrec = None
for i,line in enumerate(lines):
    rec = []
    if line[0] == 'pic':
        rec.append('Pictor A')
        picrec = i
    elif line[0] == 'cen':
        rec.append('Centaurus A')
    elif line[0] == 'vela':
        rec.append('Pup A')
    else:
        rec.append(line[0])
    rec.append(str(n.round(cat[line[0]].ra*180/n.pi,2)))
    rec.append(str(n.round(cat[line[0]].dec*180/n.pi,2)))
    spec = line[1].split(',')
    err = line[2].split(',')
    for v,e in zip(spec,err):
        rec.append(str(n.round(float(v),1)))
        rec.append(str(n.round(float(e),1)))
    if i==0:
        print "UNITS:"
        print '--'
        print 'deg'
        print 'deg'
        for i in range(len(spec)*2):
            print 'Jy'
        #column names
        print "LABLES:"
        print "Name"
        print "RAdeg"
        print "DEdeg"
        for i in range(len(spec)):
            print 'S%d'%(freqs[i])
            print 'RMS%d'%(freqs[i])
        print "EXPLANATION"
        print "Molonglo Reference Catalog (MRC) Name"
        print "J2000 RA in MRC"
        print "J2000 Dec in MRC"
        for i in range(len(spec)):
            print "beamformed flux at %d MHz"%freqs[i]
            print "rms within 10MHz channel"
    F.write(','.join(rec)+'\n')
    recs.append(rec)

#output the first 5 in a table
print """
%%note: put at top of latex file \\newcommand{\\unit}[1]{\\footnotesize \#1}
\\begin{deluxetable}{%s}
\\tablecolumns{%d}
\\tablecaption{PAPER spectra for %d MRC sources}
\\tablehead{
\\colhead{Name} &
\\colhead{\\parbox[c][3em]{2em}{Ra \\unit{deg}}} & 
\\colhead{\\parbox[c][2em]{2em}{Dec \\unit{deg}}} &"""%('l'*(len(rec)),len(rec),len(lines))
for i in range(len(spec)-1):
    print "\\colhead{\\parbox[c][2em]{2em}{S%d \\unit{Jy}}} &"%freqs[i]
    print "\\colhead{\\parbox[c][2em]{2em}{rms \\unit{Jy}}} &"
print "\\colhead{\\parbox[c][2em]{2em}{S%d \\unit{Jy}}} &"%freqs[-1]
print "\\colhead{\\parbox[c][2em]{2em}{rms \\unit{Jy}}}"

print "}"
print "\\startdata"
if not picrec is None: print '&\t'.join(recs[picrec]),'\\\\'
for i in range(10):
    print '&\t'.join(recs[i]),'\\\\'
print """
\\enddata
\\label{tab:data}
\\end{deluxetable}
"""
