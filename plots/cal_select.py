#load all the src tracks.  
from glob import glob
import numpy as n,os


def load_paper_spectra(filename):
    lines = open(filename).readlines()
    spectra ={}
    for line in lines:
        if line.startswith('#'):
            if line.startswith('#FREQ'):
                freqs = n.array(map(float,line.split('=')[1].split(','))).squeeze()
            continue
        line = line.split()
        srcname = line[0]
        fluxes = n.array(map(float,line[1].split(',')))
        errors = n.array(map(float,line[2].split(',')))
        spectra[srcname] = (freqs,fluxes,errors)
    return spectra



srctracks = glob('srctracks/*srctrack*npz')
print "found %d source tracks"%(len(srctracks))

FQ = 114
chan = 3

#load the reference calibrator first
f = n.load('srctracks/2331-416__srctrack_2455748-2455749.npz')
lst = f['lst']
i = n.argsort(lst)

trk = n.real(n.ma.masked_where(f['flags'],f['spec'])[i,FQ])
az = n.arcsin(f['x'])*180/n.pi
cal_az_count = n.ma.sum(n.abs(az)<10)


#load the uncalibrated catalog
PAPERcat  = 'psa64_pic_stripe_final.txt'
PAPER = load_paper_spectra(PAPERcat)
print "loaded %d sources from %s"%(len(PAPER),PAPERcat)
cals = []
for filename in srctracks:
    f = n.load(filename)
    srcname = os.path.basename(filename).split('__')[0]
    try:
        PAPER[srcname]
    except(KeyError):
        continue
    if srcname=='pic':continue
    i = n.argsort(f['x'])
    maxbm_possible = n.max(f['beam'][i,FQ])
    bm = n.ma.masked_where(n.logical_or(f['flags'],f['beam']==0),f['beam'])[i,FQ]
    maxbm_achieved = bm.max()
    if maxbm_achieved>=(0.75*maxbm_possible):
        cals.append(srcname)
    else:
        print srcname,maxbm_achieved
calfluxes = [n.median(PAPER[src][1]) for src in cals]
fluxorder = n.argsort(calfluxes)[::-1]
for i in fluxorder:
    flux = PAPER[cals[i]][1][chan]
    err = PAPER[cals[i]][2][chan]
    print cals[i],calfluxes[i],n.round(err/flux*100,1)
