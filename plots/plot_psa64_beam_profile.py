"""
plot psa64_beam_profile.dat with nice interpolation
"""
from matplotlib import use
use('Agg')
from pylab import *
import numpy as n
from scipy.interpolate import InterpolatedUnivariateSpline,UnivariateSpline

D = n.loadtxt('psa64_beam_profile.dat')
figure()
x = (D[:,0] - D[D[:,1].argmax(),0])*5
semilogy(x,n.abs(D[:,1]),'k')
xlabel('[arcmin]')
ylabel('psf amplitude')
subplots_adjust(bottom=0.13,left=0.13)
#Y = InterpolatedUnivariateSpline(D[:,0],D[:,1],k=5)
#x = n.linspace(D[:,0].min(),D[:,0].max())
#plot(x,Y(x),'k')
savefig('psa64_beam_profile.png')
