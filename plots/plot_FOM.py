#! /usr/bin/env python
"""
Plot the ouput of beam_src_vs_ha along with other catalog data
"""
import matplotlib
#matplotlib.use('Agg')
import os,aipy as a, numpy as n#,atpy
import optparse, sys#, scipy.optimize
#import capo as C
from pylab import *
from scipy import optimize
#import ipdb

CAT=True
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.set_description(__doc__)
#o.add_option('--mcmc_recompute',action='store_true',
#    help='Ignore saved traces and rerun mcmc fit, otherwise just move on to plotting.')
#o.add_option('--plot',action='store_true',
#    help='Do MCMC contour plots')
opts,args = o.parse_args(sys.argv[1:])

srcs,coff,cats = a.scripting.parse_srcs(opts.src,opts.cat)
cat = a.cal.get_catalog(opts.cal,srcs,catalogs=cats)
aa = a.scripting.get_null_aa()
#use PAPER for the location
aa.lat = '-30:43:17.5'
aa.lon = '21:25:41.9'
#aa.set_afreqs([0.15])
cat.compute(aa)
RAs = []
DECs = []
FOMs = []
FLUX = []
for F in args:
    lines = open(F).readlines()
    for line in lines:
        if line.startswith('#'): continue
        FOM_change = float(line.split()[-1])
        srcname = line.split()[0]
        FOMs.append(FOM_change)
        RAs.append(cat[srcname].ra)
        DECs.append(cat[srcname].dec)
        FLUX.append(cat[srcname].get_jys())
        if FOM_change>1:
            print cat[srcname],FOM_change
RAs = n.array(RAs)
DECs = n.array(DECs)
FLUX = n.array(FLUX)
figure(figsize=(7,9))
clf()
subplot(311)
plot(RAs*12/n.pi,FOMs,'.k')
xlabel('RA [h]')
ylabel('FOM improvement')
subplot(312)
plot(DECs*180/n.pi,FOMs,'.k')
xlabel('DEC [d]')
ylabel('FOM improvement')
subplot(313)
loglog(FLUX,FOMs,'.k')
xlabel('flux [Jy]')
ylabel('FOM improvement')
subplots_adjust(hspace=0.25,bottom=0.1,left=0.1)
figure()
hist(FOMs,bins=20)
xlabel('Figure of Merit increase')
show()

