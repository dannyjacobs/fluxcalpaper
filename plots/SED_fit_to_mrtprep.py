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
o.add_option('--filename',
    help='output MRT file. Default=<input>_tomrt.txt')
opts,args = o.parse_args(sys.argv[1:])

if not opts.filename is None:
    filename = opts.filename
else:
    filename = args[-1][:-4]+'_tomrt.txt'

lines = open(args[-1]).readlines()
lines = [line.split() for line in lines if not line.startswith('#')]
srclist = [l[0] for l in lines]
overlap = n.array([n.float(line[8]) for line in lines])

#if not opts.cat is None:
srclist,cutoff,catalogs = a.scripting.parse_srcs(','.join(srclist), opts.cat)
if not opts.cal is None:
    cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
else:
    cat = a.src.get_catalog(srclist, cutoff, catalogs)
o = a.scripting.get_null_aa()
cat.compute(o)

def nrow(s,n):
    return '\n'.join([s]*n)
def l2a(l,delim=','):
    #convert a delim list into a float array
    return n.array(map(float,l.split(delim)))
def p(f,d):
    "print the input as a string with n digits"
    return str(n.round(f,d))

F = open(filename,'w')
recs = []
picrec = None
total = 0
for i,line in enumerate(lines):
    rec = []
    if overlap[i]==0: continue
    if line[0] == 'pic':
        rec.append('Pictor A')
        picrec = total
    elif line[0] == 'cen':
        rec.append('Centaurus A')
    elif line[0] == 'vela':
        rec.append('Pup A')
    else:
        rec.append(line[0])

    rec.append(str(n.round(cat[line[0]].ra*180/n.pi,2)))
    rec.append(str(n.round(cat[line[0]].dec*180/n.pi,2)))
#    rec.append(line[1]) #S0
#    rec.append(line[2]) #alpha
#    rec.append(line[3]) #catalpha
#    rec.append(line[4])
    S0 = l2a(line[1])
    rec.append(p(S0[1],2))
    rec.append(p(n.abs((S0[-1]-S0[0]))/2,2))
    alpha = l2a(line[2])
    rec.append(p(alpha[1],2))
    rec.append(p(n.abs(alpha[-1] - alpha[0])/2,2))
    catS0 = l2a(line[3])
    rec.append(p(catS0[1],2))
    rec.append(p(n.abs(catS0[-1]-catS0[0])/2,2))
    catalpha = l2a(line[4])
    rec.append(p(catalpha[1],2))
    rec.append(p(n.abs(catalpha[-1] - catalpha[0])/2,2))
    
    if i==0:
        print "UNITS:"
        print '--'
        print 'deg'
        print 'deg'
        print 'Jy' #PAPER+catalog S150
        print 'Jy' #PAPER+catalog S150 76% conf
        print '--' #PAPER + catalo alpha
        print '--' #PAPER + catalo alpha 76% conf
        print 'Jy' #catalog S150 
        print 'Jy' #catalog S150 +/- 76% conf
        print '--' #catalog alpha
        print '--' #catalog alpha +/- 76% conf

        #column names
        print '\n\n'
        print "LABLES:"
        print "Name"
        print "RAdeg"
        print "DEdeg"
        print 'S150' 
        print 'e_S150' #
        print 'Sp+Index'
        print 'e_Sp+Index'
        print 'S150-prior'
        print 'e_S150-prior'
        print 'alpha-prior'
        print 'e_alpha-prior'

        print '\n\n'
        print "EXPLANATION"
        print "Molonglo Reference Catalog (MRC) Name"
        print "J2000 RA in MRC"
        print "J2000 Dec in MRC"
        print "PAPER+catalog power law fit for S150"
        print "PAPER+catalog power law fit error for S150 (76% confidence)"
        print "PAPER+catalog power law fit for Spectral Index"
        print "PAPER+catalog power law fit error for Spectral Index (76% confidence)"
        print "catalog prior power law fit for S150"
        print "catalog prior power law fit error for S150 (76% confidence)"
        print "catalog prior power law fit for Spectral Index"
        print "catalot prior power law fit error for Spectral Index (76% confidence)"

    F.write(','.join(rec)+'\n')
    recs.append(rec)
    total += 1

#output the first 5 in a table
print len(lines),total,picrec
def colhead(name,units):
    return "\\colhead{\\parbox[c][3em]{2em}{%s\\\\ \\unit{%s}}}"%(name,units)
print """
%%note: put at top of latex file \\newcommand{\\unit}[1]{\\footnotesize \#1}
\\begin{deluxetable}{%s}
\\tablecolumns{%d}
\\tablecaption{Spectral fits for %d MRC sources with and without PAPER data. 75\\%% agree well with prior
measurements. The remainder are likely contaminated by bright sources. Full table available online.}
\\tablehead{
\\multicolumn{3}{c}{ }&
\\multicolumn{4}{c}{PAPER\\tablenotemark{1} + Catalog}&
\\multicolumn{4}{c}{Catalog\\tablenotemark{2}}\\
\\colhead{Name} &
%s& 
%s &"""%('c'*len(rec),len(rec),len(recs),colhead('Ra','deg'),colhead('Dec','deg'))
print colhead('$S150$','Jy'),'&'
print colhead('$\\Delta$S','Jy'),'&'
print colhead('$\\alpha$','--'),'&'
print colhead('$\\Delta\\alpha$','--'),'&'
print colhead('$S150_p$','Jy'),'&'
print colhead('$\\Delta S_p$','Jy'),'&'
print colhead('$\\alpha_p$','--'),'&'
print colhead('$\\Delta\\alpha_p$','--')
print "}"
print "\\startdata"
if not picrec is None: print '&\t'.join(recs[picrec]),'\\\\'
for i in range(10):
    if recs[i] is None: continue
    print '&\t'.join(recs[i]),'\\\\'
print """
\\enddata
\\label{tab:fits}
\\tablenotetext{1}{Fits for the majority of sources (78\%) that agree with prior measurements are included here. 
The remainder (shown in Figure \\ref{fig:SI_contour_new}) are likely contaminated by bright sources and are not included.   Full table available online.}

\\tablenotetext{2}{MCMC fits to prior catalog data, before addition of PAPER measurements}
\\end{deluxetable}
"""
