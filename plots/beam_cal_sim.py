#! /usr/bin/env python
"""
Simulate a source track with the cal file beam, then use the bradley beam to calibrate it
What is the error?
"""


import aipy as a, numpy as n, sys, os, ephem, optparse
import matplotlib as mpl
from pylab import *
from scipy.interpolate import interp2d
def dB(x):
    return 10*n.log10(x)
def idB(x):
    return 10**(x/10.)
o = optparse.OptionParser()
o.set_usage('slice_plot.py [options]')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, src=True,pol=True)
o.add_option('--calsrc',
    help='choose a src from the list to divide out')
#o.add_option('--beamfile',
#    help='CST beam file')
opts,args = o.parse_args(sys.argv[1:])
dt = 10 /(24*60.)
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
assert(len(cat)>0)#make sure the catalog got found
if len(args)>0
    print "loading beam files ",
    print args
    files = args
    B = []
    print "freqs"
    print [f.split('_')[-1][:-4] for f in files]
    freqs = [float(f.split('_')[-1][:-4]) for f in files] #grab the frequencies in MHz
    freqs = n.array(freqs)/1e3 #convert to GHz to match with the cal file option freqs
    plotfreq = n.median(range(len(freqs))).astype(int)
    nfreqs = len(files)
    for f in files:
        B.append(n.loadtxt(f,skiprows=2))
    #reshape the datas
    B = n.array(B)
    THETA = B[:,:,0]
    THETA.shape = (nfreqs,360,181)
    PHI = B[:,:,1]
    PHI.shape = (nfreqs,360,181)
    BEAM = idB(B[:,:,2])**2 #put the beam into linear units and square
    BEAM.shape = (nfreqs,360,181)
    BEAM /= BEAM[plotfreq,0,0]
    def bradleybeam(src):
        azalt = src.get_crds('top',ncrd=2)
        theta = (n.pi/2 - azalt[1])*180/n.pi #theta is the polar coord, zero at pole
        phi = azalt[0]*180/n.pi #phi is longitude
        if opts.pol.startswith('y'): phi -= 90
        return BEAM[:,floor(phi),floor(theta)]
def rotbeam(src):
    xyz = src.get_crds('top')
    if opts.pol.startswith('x'):pol = 'y'
    else: pol = 'x'
    return aa[0].bm_response(xyz,pol=pol).squeeze()**2
def calbeam(src):
    xyz = src.get_crds('top')
    return aa[0].bm_response(xyz,pol=opts.pol[0]).squeeze()**2

plotfreq = n.median(range(len(freqs))).astype(int)
print "loading cal file"
aa = a.cal.get_aa(opts.cal,freqs)
cat.compute(aa)

#PLOT
subplot(211)
pos = {}
Abeam = {}
Bbeam = {}
figure(1)
for src in cat:
    pos[src] = []
    Abeam[src] = []
    Bbeam[src] = []
    aa.date = aa.next_rising(cat[src])
    print src
    while(True):
        t = aa.date
        aa.set_ephemtime(t+dt)
        cat.compute(aa)
        src_xyz = cat[src].get_crds('top')
        if src_xyz[2]<0:break
        src_angle = cat[src].get_crds('top',ncrd=2)
        print src_angle,bradleybeam(cat[src])[plotfreq],calbeam(cat[src])[plotfreq]
        pos[src].append(src_angle)
        Abeam[src].append(calbeam(cat[src]))
#        Bbeam[src].append(bradleybeam(cat[src]))
        Bbeam[src].append(rotbeam(cat[src]))
    pos[src] = n.array(pos[src])
    Abeam[src] = n.array(Abeam[src])
    Bbeam[src] = n.array(Bbeam[src])
    if src=='zen': color = '0.5'
    else: color = 'k'
    plot((pos[src][:,0]-n.pi)*12/n.pi,Abeam[src][:,plotfreq],color,label=src)
    plot((pos[src][:,0]-n.pi)*12/n.pi,Bbeam[src][:,plotfreq],color,linestyle=':',label=src)
xlabel('HA')
ylabel('beam amplitude')
grid()
subplot(212)
for src in cat:
    print Abeam[src].shape,Bbeam[src].shape
    bandpass = n.sum(Abeam[src] * Bbeam[src],axis=0)/n.sum(Bbeam[src]**2,axis=0)
    if not opts.calsrc is None:
        bandpass /= (n.sum(Abeam[opts.calsrc] * Bbeam[opts.calsrc],axis=0)/n.sum(Bbeam[opts.calsrc]**2,axis=0))
    if src=='zen': color='0.5'
    else: color='k'
    plot(freqs*1e3,bandpass,color,label=src)
xlabel('freqs [MHz]')
ylabel('passband amplitude')
grid()

show()
