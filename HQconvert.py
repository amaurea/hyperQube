#!/usr/bin/env python
from flipper import *
import speckMisc


def writeBbl(Bbl_T,Bbl_Cross,Bbl_Pol,dataDir,foo):
    
    numpy.savetxt(dataDir+'/BblMean_%s.dat'%foo, Bbl_T)
    numpy.savetxt(dataDir+'/BblMean_%s_Cross.dat'%foo, Bbl_Cross)
    numpy.savetxt(dataDir+'/BblMean_%s_Pol.dat'%foo, Bbl_Pol)

def compareError(l,error,error_MC,foo):
    pylab.semilogy()
    pylab.plot(l,error,label='data %s %s'%(foo,spec))
    pylab.plot(l,error_MC,label='MC %s %s'%(foo,spec))
    pylab.legend()
    pylab.show()

p = flipperDict.flipperDict()
p.read_from_file(sys.argv[1])

spec=p['spec']
dataDir='DataProduct_%s'%spec

try:
    os.makedirs(dataDir)
except:
    pass


specDir = 'spectra/'

fields=['TT','TE','EE']

for l1 in fields:
    if spec=='d5xd56':
        
        l,cl_5x5,error_5x5=numpy.loadtxt('spectra/spectrum_%s_ar1xar1_s1xs1.dat'%l1,unpack=True)
        l,cl_5x56_1,error_5x56_1=numpy.loadtxt('spectra/spectrum_%s_ar1xar1_s1xs2_a.dat'%l1,unpack=True)
        l,cl_5x56_2,error_5x56_2=numpy.loadtxt('spectra/spectrum_%s_ar1xar1_s1xs2_b.dat'%l1,unpack=True)

        Bbl_5x5_T = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_ar1xar1_s1xs1.dat')
        Bbl_5x5_Cross = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos1_ar1xar1_s1xs1.dat')
        Bbl_5x5_Pol = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos2_ar1xar1_s1xs1.dat')
        
        writeBbl(Bbl_5x5_T,Bbl_5x5_Cross,Bbl_5x5_Pol,dataDir,'d5xd5')

        Bbl_5x56_1_T = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_ar1xar1_s1xs2_a.dat')
        Bbl_5x56_1_Cross = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos1_ar1xar1_s1xs2_a.dat')
        Bbl_5x56_1_Pol = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos2_ar1xar1_s1xs2_a.dat')
        
        writeBbl(Bbl_5x56_1_T,Bbl_5x56_1_Cross,Bbl_5x56_1_Pol,dataDir,'d5xd56_p5_1')

        Bbl_5x56_2_T = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_ar1xar1_s1xs2_b.dat')
        Bbl_5x56_2_Cross = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos1_ar1xar1_s1xs2_b.dat')
        Bbl_5x56_2_Pol = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos2_ar1xar1_s1xs2_b.dat')
        
        writeBbl(Bbl_5x56_2_T,Bbl_5x56_2_Cross,Bbl_5x56_2_Pol,dataDir,'d5xd56_p5_2')

        l,cl_MC_5x5,error_MC_5x5=numpy.loadtxt('spectra_all_d5xd5/spectrum_%s_d5xd5.dat'%l1,unpack=True)
        l,cl_MC_5x56_1,error_MC_5x56_1=numpy.loadtxt('spectra_all_d5xd56_p5_1/spectrum_%s_d5xd56_p5_1.dat'%l1,unpack=True)
        l,cl_MC_5x56_2,error_MC_5x56_2=numpy.loadtxt('spectra_all_d5xd56_p5_2/spectrum_%s_d5xd56_p5_2.dat'%l1,unpack=True)
        
        compareError(l,error_5x5,error_MC_5x5,'d5xd5')
        compareError(l,error_5x56_1,error_MC_5x56_1,'d5xd56_p5_1')
        compareError(l,error_5x56_2,error_MC_5x56_2,'d5xd56_p5_2')
        
        speckMisc.writeBinnedSpectrum(l,cl_5x5,error_MC_5x5,'%s/spectrum_%s_d5xd5.dat'%(dataDir,l1))
        speckMisc.writeBinnedSpectrum(l,cl_5x56_1,error_MC_5x56_1,'%s/spectrum_%s_d5xd56_p5_1.dat'%(dataDir,l1))
        speckMisc.writeBinnedSpectrum(l,cl_5x56_2,error_MC_5x56_2,'%s/spectrum_%s_d5xd56_p5_2.dat'%(dataDir,l1))

    if spec=='d6xd56':

        l,cl_6x6,error_6x6=numpy.loadtxt('spectra/spectrum_%s_ar1xar1_s1xs1.dat'%l1,unpack=True)
        l,cl_6x56_1,error_6x56_1=numpy.loadtxt('spectra/spectrum_%s_ar1xar1_s1xs2_a.dat'%l1,unpack=True)
        l,cl_6x56_2,error_6x56_2=numpy.loadtxt('spectra/spectrum_%s_ar1xar1_s1xs2_b.dat'%l1,unpack=True)

        Bbl_6x6_T = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_ar1xar1_s1xs1.dat')
        Bbl_6x6_Cross = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos1_ar1xar1_s1xs1.dat')
        Bbl_6x6_Pol = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos2_ar1xar1_s1xs1.dat')

        writeBbl(Bbl_6x6_T,Bbl_6x6_Cross,Bbl_6x6_Pol,dataDir,'d6xd6')

        Bbl_6x56_1_T = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_ar1xar1_s1xs2_a.dat')
        Bbl_6x56_1_Cross = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos1_ar1xar1_s1xs2_a.dat')
        Bbl_6x56_1_Pol = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos2_ar1xar1_s1xs2_a.dat')

        writeBbl(Bbl_6x56_1_T,Bbl_6x56_1_Cross,Bbl_6x56_1_Pol,dataDir,'d6xd56_p6_1')

        Bbl_6x56_2_T = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_ar1xar1_s1xs2_b.dat')
        Bbl_6x56_2_Cross = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos1_ar1xar1_s1xs2_b.dat')
        Bbl_6x56_2_Pol = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos2_ar1xar1_s1xs2_b.dat')

        writeBbl(Bbl_6x56_2_T,Bbl_6x56_2_Cross,Bbl_6x56_2_Pol,dataDir,'d6xd56_p6_2')

        l,cl_MC_6x6,error_MC_6x6=numpy.loadtxt('spectra_all_d6xd6/spectrum_%s_d6xd6.dat'%l1,unpack=True)
        l,cl_MC_6x56_1,error_MC_6x56_1=numpy.loadtxt('spectra_all_d6xd56_p6_1/spectrum_%s_d6xd56_p6_1.dat'%l1,unpack=True)
        l,cl_MC_6x56_2,error_MC_6x56_2=numpy.loadtxt('spectra_all_d6xd56_p6_2/spectrum_%s_d6xd56_p6_2.dat'%l1,unpack=True)

        compareError(l,error_6x6,error_MC_6x6,'d6xd6')
        compareError(l,error_6x56_1,error_MC_6x56_1,'d6xd56_p6_1')
        compareError(l,error_6x56_2,error_MC_6x56_2,'d6xd56_p6_2')

        speckMisc.writeBinnedSpectrum(l,cl_6x6,error_MC_6x6,'%s/spectrum_%s_d6xd6.dat'%(dataDir,l1))
        speckMisc.writeBinnedSpectrum(l,cl_6x56_1,error_MC_6x56_1,'%s/spectrum_%s_d6xd56_p6_1.dat'%(dataDir,l1))
        speckMisc.writeBinnedSpectrum(l,cl_6x56_2,error_MC_6x56_2,'%s/spectrum_%s_d6xd56_p6_2.dat'%(dataDir,l1))


    if spec=='d56_1xd56_1':
        l,cl_56_1x56_1,error_56_1x56_1=numpy.loadtxt('spectra/spectrum_%s_ar1xar1_s2xs2.dat'%l1,unpack=True)
        Bbl_56_1x56_1_T = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_ar1xar1_s2xs2.dat')
        Bbl_56_1x56_1_Cross = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos1_ar1xar1_s2xs2.dat')
        Bbl_56_1x56_1_Pol = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos2_ar1xar1_s2xs2.dat')

        writeBbl(Bbl_56_1x56_1_T,Bbl_56_1x56_1_Cross,Bbl_56_1x56_1_Pol,dataDir,'d56_1xd56_1')
        
        j,cl_MC_56_1x56_1,error_MC_56_1x56_1=numpy.loadtxt('spectra_all_d56_1xd56_1/spectrum_%s_d56_1xd56_1.dat'%l1,unpack=True)
        speckMisc.writeBinnedSpectrum(l,cl_56_1x56_1,error_MC_56_1x56_1,'%s/spectrum_%s_d56_1xd56_1.dat'%(dataDir,l1))

    if spec=='d56_1xd56_2':
        l,cl_56_1x56_2,error_56_1x56_2=numpy.loadtxt('spectra/spectrum_%s_ar1xar2_s2xs2.dat'%l1,unpack=True)
        Bbl_56_1x56_2_T = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_ar1xar2_s2xs2.dat')
        Bbl_56_1x56_2_Cross = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos1_ar1xar2_s2xs2.dat')
        Bbl_56_1x56_2_Pol = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos2_ar1xar2_s2xs2.dat')
    
        writeBbl(Bbl_56_1x56_2_T,Bbl_56_1x56_2_Cross,Bbl_56_1x56_2_Pol,dataDir,'d56_1xd56_2')
    
        j,cl_MC_56_1x56_2,error_MC_56_1x56_2=numpy.loadtxt('spectra_all_d56_1xd56_2/spectrum_%s_d56_1xd56_2.dat'%l1,unpack=True)
        speckMisc.writeBinnedSpectrum(l,cl_56_1x56_2,error_MC_56_1x56_2,'%s/spectrum_%s_d56_1xd56_2.dat'%(dataDir,l1))

    if spec=='d56_2xd56_2':
        l,cl_56_2x56_2,error_56_2x56_2=numpy.loadtxt('spectra/spectrum_%s_ar2xar2_s2xs2.dat'%l1,unpack=True)
        Bbl_56_2x56_2_T = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_ar2xar2_s2xs2.dat')
        Bbl_56_2x56_2_Cross = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos1_ar2xar2_s2xs2.dat')
        Bbl_56_2x56_2_Pol = numpy.loadtxt('mcm/Bbl.mcm_mask10amin_cos2_ar2xar2_s2xs2.dat')
    
        writeBbl(Bbl_56_2x56_2_T,Bbl_56_2x56_2_Cross,Bbl_56_2x56_2_Pol,dataDir,'d56_2xd56_2')
    
        j,cl_MC_56_2x56_2,error_MC_56_2x56_2=numpy.loadtxt('spectra_all_d56_2xd56_2/spectrum_%s_d56_2xd56_2.dat'%l1,unpack=True)
        speckMisc.writeBinnedSpectrum(l,cl_56_2x56_2,error_MC_56_2x56_2,'%s/spectrum_%s_d56_2xd56_2.dat'%(dataDir,l1))

