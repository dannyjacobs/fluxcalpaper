"""
python script to sort out output from fit_SED.py
"""
import numpy as n,sys
from pylab import *

lines = open(sys.argv[-1]).readlines()
lines = [line for line in lines if not line.startswith('#')]
print [len(line.split()) for line in lines]
improvement = n.array([n.float(line.split()[9]) for line in lines]) # (FOM change)*(fraction of PAPER overlapping previos)
overlap = n.array([n.float(line.split()[8]) for line in lines])
coastline = n.array([float(line.split()[10]) for line in lines]) # number of pixels on the 76% contour
poor = []
change = []
j = 0
coastmax = n.mean(coastline)+n.std(coastline)*3
disagree = []
imps = []
goodlines = []
shit = []
for i in n.argsort(improvement)[::-1]:
    if coastline[i]>coastmax: 
        #what else could there be?  
        print '%Bad fit: ',lines[i].strip()
        print '%\t\t improvement = ',improvement[i],'\t coastline = ',coastline[i]
        shit.append(line)
        continue #throw out bad fits with bubbly contours (long "coastline")  As of 5 April 2013 there are none!
    if improvement[i] >0 and overlap[i]>0.001:   
        goodlines.append(lines[i]) 
        #good safe fits have some overlap, positive improvement and not bubbly contours (coastline) #on April 5 there
        #is just one of these! I'm probably going to just include it in the good fits.
    elif improvement[i]<0 and overlap[i]>0.001:
        poor.append(lines[i])
        #actually negative improvement 
    else:
        #these should just be sources that disagree completely
        disagree.append(lines[i])


print  "\\begin{figure*}[htbp]"
print  "\\begin{center}"
print "%ok fit, but doesn't agree with past data, or actually increases model uncertainty"
print "%contour length limit = ",coastmax
for line in n.sort(poor):
    print '%Poor improvement: ',line.strip()
    print "\includegraphics[width=2in]{plots/%s_SI_MCMC.png} %%%s"%\
    (line.split()[0],line.split()[9]) #
    

print """\\end{center}
\\caption{fits of the worst sources. Contours as described in Figure \\ref{fig:SI_contour_1}. The fit either disagrees
with past data or does not offer any improvement (improvement index $\le 0$). 0008-421 is consistent with evidence for obsorption, the rest are near
the galactic plane and are likely dominated by sidelobes.
}\\label{fig:SI_contour_bad}
\\end{figure*}"""

print "\\begin{figure*}[htbp]"
print "\\begin{center}"
print "%zero overlap"
for line in n.sort(disagree):
    print "\\includegraphics[width=2in]{plots/%s_SI_MCMC.png} %%%s"%(line.split()[0],line.split()[9])
print """\\end{center}
\\caption{These sources offer results at odds with previous measurements.}
\\label{fig:SI_contour_new}
\\end{figure*}
"""


#for j,line in enumerate(goodlines):
#    print line.strip()
perfig=12
print "%the best sources. High quality fit that improves on previous knowledge" 
print  "\\begin{figure*}[htbp]"
print  "\\begin{center}"
improvement = n.array([n.float(line.split()[9]) for line in n.sort(goodlines)])
for j,line in enumerate(n.sort(goodlines)):
    #print lines[i].strip()
    if line.split()[0]=='pic': 
        picimprovement=improvement[j]
        pic = line
        continue
    if not j%perfig and j!=0:
        print """\\end{center}
\\caption{Spectral model contours as described in Figure \\ref{fig:pic_spectrum}. Sources marked with a
* were used to assess calibration error.
}\\label{fig:SI_contour_%d}
\\end{figure*}
\\clearpage
\\begin{figure*}[htbp]2\\begin{center}"""%n.ceil(j/perfig)
    print "\includegraphics[width=2in]{plots/%s_SI_MCMC.png} %%%s"%(line.split()[0],line.split()[9])



print """\\end{center}
\\caption{fits of the next 16 sources, as described in Figure \\ref{fig:SI_contour_1}. Sources marked with a
* were used to assess calibration error.
}\\label{fig:SI_contour_%d}
\\end{figure*}
"""%(n.ceil(j/perfig)+1)
#cnts,bins = histogram(improvement,bins=n.logspace(n.log10(.1),n.log10(4),num=10))
#semilogx(bins[1:],cnts,'k')
#vlines(0.4,cnts.min(),cnts.max())
improvement = n.array([n.float(line.split()[9]) for line in lines])
print "%%SUMMARY"
print "%% wtf shit:",len(shit)
print "%%>0 improvement:",len(goodlines)
print "%%<0 improvement:",len(poor)
print "%%=0 improvement:",len(disagree)
print "%%total %d"%len(lines)
print "%%median improvement: ",n.median(improvement)
print "%%max improvement:",n.max(improvement)
print "%%min improvement:",n.min(improvement)
print '%%sources "better" than pictor:',n.sum(improvement>picimprovement),'/',n.sum(improvement>0)
print '%%',pic
for i,line in enumerate(n.sort(goodlines)):
    if improvement[i]>(picimprovement):
        print '%%',line.strip()
show()
