#! /usr/bin/env python
"""
Simulate a source track with the cal file beam, then use the bradley beam to calibrate it
What is the error?
"""
#./beam_cal_sim.py -C psa64_CSTbeam -s pic,2331-416,zen -pyy --cat=southern_sky_v3 --calsrc=2331-416 /Users/danny/Work/radio_astronomy/2010_beam/sdipole_05e_eg_ffx_*.txt

import aipy as a, numpy as n, sys, os, ephem, optparse
import matplotlib as mpl
from pylab import *
from glob import glob
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
opts,args = o.parse_args(sys.argv[1:])
dt = 10 /(24*60.)
srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
assert(len(cat)>0)#make sure the catalog got found
print "loading CST beam files ",
#files = [f for f in args if f.endswith('.txt')]
files  = glob('sdipole*ffx*150.txt')
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
def rotbradleybeam(src):
    azalt = src.get_crds('top',ncrd=2)
    theta = (n.pi/2 - azalt[1])*180/n.pi #theta is the polar coord, zero at pole
    phi = azalt[0]*180/n.pi #phi is longitude
    if opts.pol.startswith('y'): phi -= 90
    phi -= 90
    return BEAM[:,floor(phi),floor(theta)]
    
print "Loading fits beam"
jcp_beam = a.map.Map()
jcp_beam.from_fits('test19b_f6_sm.fits')

def rotbeam(src):
    xyz = src.get_crds('top')
    if opts.pol.startswith('x'):pol = 'y'
    else: pol = 'x'
    return aa[0].bm_response(xyz,pol=pol).squeeze()**2
def calbeam(src):
    xyz = src.get_crds('top')
    bm = aa[0].bm_response(xyz,pol=opts.pol[0]).squeeze()**2
    if type(bm)!=n.ndarray: bm = n.array([bm])
    return bm
def ModelBeam(src):
    return calbeam(src)
def TrueBeam(src):
#    return rotbradleybeam(src)
#    return calbeam(src)
    return jcp_beam[tuple(src.get_crds('top'))]**2
nfreqs = len(freqs)
plotfreq = int(nfreqs/2)
print "loading cal file"
aa = a.cal.get_aa(opts.cal,freqs)
cat.compute(aa)

#PLOT
subplot(211)
pos = {}
Abeam = {}
Bbeam = {}

print "Beam model amplitudes at meridian"
for src in cat:
    aa.date = aa.next_transit(cat[src])
    cat.compute(aa)
    print "beamfreq=",aa.get_afreqs()[0]
    print src,cat[src].get_crds('top',ncrd=2),
    print ModelBeam(cat[src])[0],TrueBeam(cat[src])[0]


#scan through a range of decs and simulate the weight factor g using two different beams
decs = n.arange(-59,-30,1)
srclist = ','.join(['0_%d'%dec for dec in decs])
srcs,coff,catalogs = a.scripting.parse_srcs(srclist,'')
cat = a.src.get_catalog(srcs,catalogs)
cat.compute(aa)
g = n.zeros((nfreqs,len(cat)))
w = n.zeros((nfreqs,len(cat)))
A = {}
B = {}
A2 = {}
W = {}
crds = {}
times = {}
for i,src in enumerate(n.sort(cat.keys())[::-1]):
    print src
    aa.date = aa.next_rising(cat[src])
    A[src] = []
    B[src] = []
    A2[src] = []
    crds[src] = []
    W[src] = n.zeros(nfreqs)
    times[src] = []
    while True:
        t = aa.date
        aa.set_ephemtime(t+dt)
        cat.compute(aa)        
        src_xyz = cat[src].get_crds('top')
        src_ang = cat[src].get_crds('top',ncrd=2)
        if src_xyz[2]<0:break
        if n.abs(src_ang[1]-n.pi)<(10*n.pi/180): break
        print src_ang
        crds[src].append(src_xyz)
        g[:,i] += ModelBeam(cat[src]) * TrueBeam(cat[src])
        w[:,i] += ModelBeam(cat[src]) * ModelBeam(cat[src])
        A[src].append(ModelBeam(cat[src]))
        B[src].append(TrueBeam(cat[src]))
        A2[src].append(bradleybeam(cat[src]))
        W[src] += ModelBeam(cat[src]) * ModelBeam(cat[src])
        times[src].append(t)
    A[src] = n.array(A[src])
    B[src] = n.array(B[src])
    A2[src] = n.array(A2[src])

g /= w
figure()
clf()
caldec = 15
print "calibrating to declination %d"%decs[caldec]
for src in ['0_-40','0_-50']:
    t = n.array(times[src])
    t -= n.median(t)
    t *= 24
    B_A = n.array(A[src])
    B_M = n.array(B[src])
    subplot(211)
    plot(t,B_A[:,plotfreq])
    plot(t,B_M[:,plotfreq],ls=':')
    subplot(212)
    plot(t,B_A[:,plotfreq]*B_M[:,plotfreq])
subplots_adjust(bottom=0.1,left=0.1)
xlabel('time [hrs]')
ylabel('$B_M(\\theta,\\phi) B_A(\\theta,\\phi)$')
print "Model weighted, simulated beam"
savefig('model_weighted_beam_sim_track.png')
print "\t model_weighted_beam_sim_track.png"
figure()
#subplot(211)
for i in range(len(decs)):
    plot(freqs*1e3,g[:,i]/g[:,caldec],'0.5')
grid()
xlabel('freq [MHz]')
ylabel('g($\\delta$)')

#compare the two beam models
g2 = n.array([n.sum(A2['0_%d'%dec]*B['0_%d'%dec])/n.sum(A2['0_%d'%dec]**2) for dec in decs])
print g2.shape
figure()
#subplot(212)
plot(decs,g[plotfreq,:]/g[plotfreq,caldec],'k')
plot(decs,g2/g2[caldec],'k')
subplots_adjust(bottom=0.1,left=0.1)
xlabel('$\\delta$  [deg]')
ylabel('g($\\delta$)')
grid()
savefig('beam_gain_sim.png')
print "Model-weighted beam integration sim. passband and declination dependence"
print "\t beam_gain_sim.png"

#show()
