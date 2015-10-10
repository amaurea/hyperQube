#!/usr/bin/env python
from flipper import *
import speckMisc




p = flipperDict.flipperDict()
p.read_from_file(sys.argv[1])
dir=sys.argv[2]

os.system('cp -r %s/* .'%dir)

arrayTags = p['arrayTags']
seasonTags = p['seasonTags']

array=arrayTags[0]
season=seasonTags[0]


patchDir = 'patches'
l = os.listdir(patchDir)
nDivs = 0
nPatches = 0

for il in l:
    if 'all' in il:
        continue
    if 'T_map_%s_%s'%(arrayTags[0],seasonTags[0]) in il:
        nDivs += 1
    if 'T_map_%s_%s_0'%(arrayTags[0],seasonTags[0]) in il and '_0' in il[-2:]:
        nPatches += 1




n=100
phiMax=2
ang=numpy.linspace(-phiMax,phiMax,n)
chi2=numpy.zeros(n)

count=0
for phi in ang:
    phi*=numpy.pi/(180)
    
    for i in xrange(nDivs):
            
        T = liteMap.liteMapFromFits('%s/%s/T_map_%s_%s_%d'%(dir,patchDir,array,season,i))
        Q = liteMap.liteMapFromFits('%s/%s/Q_map_%s_%s_%d'%(dir,patchDir,array,season,i))
        U = liteMap.liteMapFromFits('%s/%s/U_map_%s_%s_%d'%(dir,patchDir,array,season,i))
        
        Q_rot=Q.copy()
        U_rot=U.copy()
        
        Q_rot.data=Q.data*numpy.cos(2*phi)+U.data*numpy.sin(2*phi)
        U_rot.data=-Q.data*numpy.sin(2*phi)+U.data*numpy.cos(2*phi)

        Q_rot.writeFits(patchDir+os.path.sep+'Q_map_%s_%s_%d'%(array,season,i),overWrite=True)
        U_rot.writeFits(patchDir+os.path.sep+'U_map_%s_%s_%d'%(array,season,i),overWrite=True)

    os.system('HQcompileSpectra.py global.dict')
    os.system('HQcomputeAnalyticCovariance.py global.dict')

    print 'only seasonTags[0] arrayTags[0] at the moment'
    l,cl_EB,error_EB=numpy.loadtxt('spectra/spectrum_EB_%sx%s_%sx%s.dat'%(array,array,season,season),unpack=True)


    speckMisc.writeBinnedSpectrum(l,cl_EB,error_EB,'spectrum_EB_%sx%s_%sx%s_%d.dat'%(array,array,season,season,count))
    chi2[count]=numpy.mean(cl_EB**2/error_EB**2)

    print phi,chi2[count]
    count+=1

pylab.plot(ang,chi2)
pylab.show()