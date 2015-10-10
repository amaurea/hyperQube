#!/usr/bin/env python


from flipper import *
#from scipy.integrate import 
from scipy.interpolate import splev,splrep
from flipper import *
from numpy.fft import fftshift,fftfreq,fft2,ifft2, ifftshift
from scipy import interpolate
from scipy import *
import scipy
import os
import random
import sys
import pickle
import time
import fortran



def makeTemplate(m, wl, ell, maxEll, outputFile = None):
    """
    Yanked from Toby's csFilter
    For a given map (m) return a 2D k-space template from a 1D specification wl
    ell = 2pi * i / deltaX
    (m is not overwritten)
    """
    
    ell = numpy.array(ell)
    wl  = numpy.array(wl)
    
    p2d = fftTools.powerFromLiteMap(m)
    p2d.powerMap[:] = 0.
    # print "max_lx, max_ly", p2d.lx.max(), p2d.ly.max()
    # print "m_dx, m_dy", m.pixScaleX, m.pixScaleY
    # print "m_nx, m_ny", m.Nx, m.Ny
    l_f = numpy.floor(p2d.modLMap)
    l_c = numpy.ceil(p2d.modLMap)
    
    for i in xrange(numpy.shape(p2d.powerMap)[0]):
        for j in xrange(numpy.shape(p2d.powerMap)[1]):
            if l_f[i,j] > maxEll or l_c[i,j] > maxEll:
                continue
            w_lo = wl[l_f[i,j]]
            w_hi = wl[l_c[i,j]]
            trueL = p2d.modLMap[i,j]
            w = (w_hi-w_lo)*(trueL - l_f[i,j]) + w_lo
            p2d.powerMap[i,j] = w
            
    if outputFile != None:
        p2d.writeFits(outputFile, overWrite = True)
    return p2d


def trimShiftKMap(kMap,elTrim,Lx,Ly,lx,ly,type):
    """
    @brief Trims a 2-D powerMap at a certain Lx, Ly, and returns the trimmed kMap object
    trimmed at elTrim and shifted by Lx and Ly
    """
    # assert(elTrim>0.)
    
    # assert((lx[0] != 0.) and (ly[0] != 0.))
    #kM = kMap.copy()
    idx = numpy.where((lx < elTrim+Lx) & (lx > -elTrim+Lx))
    idy = numpy.where((ly < elTrim+Ly) & (ly > -elTrim+Ly))
    
    # print "numpy.where takes %f secs"%(time.time()-ta)
    
    trimA = kMap[idy[0],:]
    trimB = trimA[:,idx[0]]
    kTrimmed = trimB/(numpy.sqrt(float(Lx)**2+float(Ly)**2)+2)**(type*2)

    #print t1,t2,t3
    
    #del trimA,trimB
    # print "plus trim  takes %f secs"%(time.time()-ta)
    
    # trimB.flags
    kTrimmed=numpy.ascontiguousarray(trimB)
    # trimB.flags
    
    return kTrimmed


def generateMCM(window1,window2,binLo,binHi, trimAtL=10000,mbbFilename = None,transfer = None,  binningWeightMap = None, type=0,std=0,pixWin=True,arrayTag=None,spTag=None):
    """
    window: data window
    
    """
    
    print "TYPE:",type
    
    powerOfL=0
    
    powerMaskHolder_11 = fftTools.powerFromLiteMap(window1,window1)
    
    pixW_temp=powerMaskHolder_11.copy()
    pixW = powerMaskHolder_11.pixelWindow()
    pixW_temp.powerMap[:]=pixW
    pixW_temp= pixW_temp.trimAtL(2*trimAtL+500)

    powerMaskHolder_11= powerMaskHolder_11.trimAtL(2*trimAtL+500)
    powerMaskHolder_12 = fftTools.powerFromLiteMap(window1,window2)
    powerMaskHolder_12 = powerMaskHolder_12.trimAtL(2*trimAtL+500)

    powerMaskHolder_22 = fftTools.powerFromLiteMap(window2,window2)
    powerMaskHolder_22= powerMaskHolder_22.trimAtL(2*trimAtL+500)

    
    
    powerMask_11 = powerMaskHolder_11.powerMap.copy()
    powerMaskShifted_11 = fftshift(powerMask_11)
    
    powerMask_12  = powerMaskHolder_12.powerMap.copy()
    powerMaskShifted_12  = fftshift(powerMask_12 )
    
    powerMask_22 = powerMaskHolder_22.powerMap.copy()
    powerMaskShifted_22 = fftshift(powerMask_22)
    
    phlx = fftshift(powerMaskHolder_11.lx)
    phly = fftshift(powerMaskHolder_11.ly)
    
    
    if transfer != None:
        
        if len(transfer)==4:
            ell,f_ell_Pol,f_ell_Cross,f_ell_TT=transfer
        else:
            ell,f_ell_Pol= transfer
            f_ell_TT = f_ell_Pol.copy()
            f_ell_Cross= f_ell_Pol.copy()
        
        t_TT = makeTemplate( window1.copy(), f_ell_TT, ell, trimAtL )
        t_Cross = makeTemplate( window1.copy(), f_ell_Cross, ell, trimAtL )
        t_Pol = makeTemplate( window1.copy(), f_ell_Pol, ell, trimAtL )

        if pixWin!=False:
            print 'Apply Pix Window Function'
            t_TT.powerMap *= pixW
            t_Cross.powerMap *= pixW
            t_Pol.powerMap *= pixW
        else:
            print 'No pix Win applied'
        
    
        transferTrim_TT = t_TT.trimAtL(trimAtL)
        transferTrim_Cross = t_Cross.trimAtL(trimAtL)
        transferTrim_Pol = t_Pol.trimAtL(trimAtL)

        
    if binningWeightMap ==None:
        binningWeightMap = powerMask_11.copy()*0.+1.0
    else:
        assert(binningWeightMap.shape == powerMask_11.shape)
        
        
    powerMaskHolderTrim = powerMaskHolder_11.trimAtL(trimAtL)
    ly=powerMaskHolderTrim.ly
    lx=powerMaskHolderTrim.lx
    iy, ix = numpy.mgrid[0:powerMaskHolderTrim.Ny,0:powerMaskHolderTrim.Nx]
    
    powerMaskTrim = powerMaskHolderTrim.powerMap.copy()
    print "Trimmed kMap dims",powerMaskHolderTrim.Ny,powerMaskHolderTrim.Nx
    assert(binHi.max() <= trimAtL)
    
    pMMaps = []
    pMMaps_cos = []
    pMMaps_sin = []
    pMMaps_cos2 = []
    pMMaps_sin2 = []
    pMMaps_cossin = []

    #bins = [[0.,200.],[200.,400.],[400.,700.],[700.,1200.],[1200,1800],[1800,2500],[2500,3500]]

    mArray = numpy.zeros(shape=(len(binLo),len(binLo)))
    mArray_cos = mArray.copy()
    mArray_sin = mArray.copy()
    mArray_cos2 = mArray.copy()
    mArray_sin2 = mArray.copy()
    mArray_cossin = mArray.copy()

    Bbl = numpy.zeros(shape=(len(binLo),numpy.int(trimAtL)))
    Bbl_cos = Bbl.copy()
    Bbl_sin = Bbl.copy()
    Bbl_cos2 = Bbl.copy()
    Bbl_sin2 = Bbl.copy()
    Bbl_cossin = Bbl.copy()

    modIntLMap = numpy.array(powerMaskHolder_11.modLMap + 0.5,dtype='int64')
    t00 = time.time()
    cumTerms = 0
  
  	
    ang=numpy.fft.fftshift(numpy.arctan2(ly[iy],lx[ix]))
    cos_array=numpy.cos(2*ang)
    sin_array=numpy.sin(2*ang)


    for ibin in xrange(len(binLo)):
        t0 = time.time()
        
        location = numpy.where((modIntLMap >= binLo[ibin]) & (modIntLMap <= binHi[ibin]))
        
        #nMax=2000
        #print len(location[0])
        #if len(location[0])>nMax:
            
            #    id=numpy.random.randint(0,len(location[0]),nMax)
            #id=list(set(id))
            #location=[location[0][id],location[1][id]]

#print len(location[0])
        
        print "where on modIntLMap took %f secs"%(time.time()-t0)
        binMap = powerMask_11.copy()*0.
        binMap[location] = binningWeightMap[location]
        sumBin = binMap.sum()
        binMap[location] *= powerMaskHolder_11.modLMap[location]**powerOfL
        #print binningWeightMap[location]
        assert(sumBin>0.)
        
        binMap0 = (trimShiftKMap(powerMaskShifted_11,trimAtL,0,0,phlx,phly,type))*0.
        binMap_cos = binMap0.copy()
        binMap_sin = binMap0.copy()
        binMap_cos2 = binMap0.copy()
        binMap_sin2 = binMap0.copy()
        binMap_cossin = binMap0.copy()

        
        t0 = time.time()
        deltaT = 0.
        cumTerms += len(location[0])
 # computeBinMap(location[0],location[1],powerMaskHolder_11.ly,powerMaskHolder_11.lx,powerMaskShifted_11,powerMaskShifted_12,powerMaskShifted_22,cos_array,sin_array,binMap,binMap0, binMap_cos,binMap_sin,binMap_cos2,binMap_sin2,binMap_cossin,phlx,phly,type)

        fortran.fortran(location[0],location[1],powerMaskHolder_11.ly,powerMaskHolder_11.lx,powerMaskShifted_11.T,powerMaskShifted_12.T,powerMaskShifted_22.T,cos_array.T,sin_array.T,binMap.T,binMap0.T, binMap_cos.T,binMap_sin.T,binMap_cos2.T,binMap_sin2.T,binMap_cossin.T,phlx,phly,type,trimAtL)


#print "Total time in trimshift calls %f secs"%deltaT
        print "Time to do sum  of %d terms %f secs"%(len(location[0]),time.time() - t0)
        print "Best estimate of time per step %f"%((time.time() - t00)/cumTerms)
        
        binMap0 = ifftshift(binMap0)
        binMap_cos = ifftshift(binMap_cos)
        binMap_sin = ifftshift(binMap_sin)
        binMap_cos2 = ifftshift(binMap_cos2)
        binMap_sin2 = ifftshift(binMap_sin2)
        binMap_cossin = ifftshift(binMap_cossin)
        
        
        if transfer != None:
            pMMaps.append(binMap0[:]*transferTrim_TT.powerMap[:]/sumBin)
            pMMaps_cos.append(binMap_cos[:]*transferTrim_Cross.powerMap[:]/sumBin)
            pMMaps_sin.append(binMap_sin[:]*transferTrim_Cross.powerMap[:]/sumBin)
            pMMaps_cos2.append(binMap_cos2[:]*transferTrim_Pol.powerMap[:]/sumBin)
            pMMaps_sin2.append(binMap_sin2[:]*transferTrim_Pol.powerMap[:]/sumBin)
            pMMaps_cossin.append(binMap_cossin[:]*transferTrim_Pol.powerMap[:]/sumBin)
        else:
            pMMaps.append(binMap0/sumBin)
            pMMaps_cos.append(binMap_cos/sumBin)
            pMMaps_sin.append(binMap_sin/sumBin)
            pMMaps_cos2.append(binMap_cos2/sumBin)
            pMMaps_sin2.append(binMap_sin2/sumBin)
            pMMaps_cossin.append(binMap_cossin/sumBin)
        
        print 'Bin %d (%f,%f) of %d bins took %f secs to compute'\
              %(ibin, binLo[ibin],binHi[ibin],len(binLo),(time.time() - t0))
    
    print "PM maps done in %f secs"%(time.time() - t00)
        
    larray = numpy.arange(2,numpy.int(trimAtL)+2)
    deltaLx = numpy.abs(powerMaskHolderTrim.modLMap[0,1] - powerMaskHolderTrim.modLMap[0,0])
    deltaLy = numpy.abs(powerMaskHolderTrim.modLMap[1,0] - powerMaskHolderTrim.modLMap[0,0])
    delta = numpy.min([deltaLx,deltaLy])/2.0
    gaussFactor = 1./numpy.sqrt(2.*numpy.pi*delta**2.)

    modLMapInt = numpy.array(powerMaskHolderTrim.modLMap+0.5, dtype='int64')

    t0 = time.time()
    
    for j in xrange(len(binHi)):
        location = numpy.where((modLMapInt >= binLo[j]) &\
                               (modLMapInt <= binHi[j]))
        binMapTrim = powerMaskTrim.copy()*0.
        binMapTrim[location] = 1.
        binMapTrim[location] *= numpy.nan_to_num(1./(powerMaskHolderTrim.modLMap[location])**powerOfL)
        for i in xrange(len(pMMaps)):
            newMap = pMMaps[i].copy()
            newMap_cos =pMMaps_cos[i].copy()
            newMap_sin =pMMaps_sin[i].copy()
            newMap_cos2 =pMMaps_cos2[i].copy()
            newMap_sin2 =pMMaps_sin2[i].copy()
            newMap_cossin =pMMaps_cossin[i].copy()
            
            result = (newMap*binMapTrim).sum()
            result_cos = (newMap_cos*binMapTrim).sum()
            result_sin = (newMap_sin*binMapTrim).sum()
            result_cos2 = (newMap_cos2*binMapTrim).sum()
            result_sin2 = (newMap_sin2*binMapTrim).sum()
            result_cossin = (newMap_cossin*binMapTrim).sum()
            
            
            mArray[i,j] =result
            mArray_cos[i,j] = result_cos
            mArray_sin[i,j] = result_sin
            mArray_cos2[i,j] = result_cos2
            mArray_sin2[i,j] = result_sin2
            mArray_cossin[i,j] = result_cossin

    print "MArray done in %f secs"%(time.time()-t0)
        
    modMap = numpy.ravel(powerMaskHolderTrim.modLMap)    
#modMap=  numpy.floor(modMap)
    deltaLx = numpy.abs(powerMaskHolderTrim.lx[1]-powerMaskHolderTrim.lx[0])
    deltaLy = numpy.abs(powerMaskHolderTrim.ly[1]-powerMaskHolderTrim.ly[0])
    

    t0 = time.time()
    for k in xrange(len(larray)):
        
        gauss = numpy.exp(-(larray-larray[k])**2./(2.*delta**2.))
        sum = gauss.sum()

        gauss = numpy.exp(-(modMap-larray[k])**2./(2.*delta**2.))
        gauss /= sum
        binMapTrim = numpy.reshape(gauss,[powerMaskHolderTrim.Ny,powerMaskHolderTrim.Nx])
                
        binMapTrim *= powerMaskHolderTrim.modLMap**powerOfL
        for i in xrange(len(pMMaps)):
            result = (pMMaps[i]*binMapTrim).sum()
            result_cos = (pMMaps_cos[i]*binMapTrim).sum()
            result_sin = (pMMaps_sin[i]*binMapTrim).sum()
            result_cos2 = (pMMaps_cos2[i]*binMapTrim).sum()
            result_sin2 = (pMMaps_sin2[i]*binMapTrim).sum()
            result_cossin = (pMMaps_cossin[i]*binMapTrim).sum()
            
            
            Bbl[i,k] = result
            Bbl_cos[i,k] = result_cos
            Bbl_sin[i,k] = result_sin
            Bbl_cos2[i,k] = result_cos2
            Bbl_sin2[i,k] = result_sin2
            Bbl_cossin[i,k] = result_cossin

    print "Bbl done in %f secs"%(time.time()-t0)
    
    Bbl = numpy.dot(scipy.linalg.inv(mArray), Bbl)
    Bbl_cos = numpy.dot(scipy.linalg.inv(mArray_cos), Bbl_cos)
    Bbl_sin = numpy.dot(scipy.linalg.inv(mArray_sin), Bbl_sin)
    Bbl_cos2 = numpy.dot(scipy.linalg.inv(mArray_cos2), Bbl_cos2)
    Bbl_sin2 = numpy.dot(scipy.linalg.inv(mArray_sin2), Bbl_sin2)
    Bbl_cossin = numpy.dot(scipy.linalg.inv(mArray_cossin), Bbl_cossin)

    if mbbFilename != None:
        mbbDir = (mbbFilename.split("/"))[0]+"/"
        mbbFilename = (mbbFilename.split("/"))[-1]
        numpy.savetxt('%s%s_%s_%s.dat'%(mbbDir,mbbFilename,arrayTag,spTag),mArray)
        numpy.savetxt('%s%s_%s_%s.dat'%(mbbDir,mbbFilename+'_cos1',arrayTag,spTag),mArray_cos)
        numpy.savetxt('%s%s_%s_%s.dat'%(mbbDir,mbbFilename+'_sin1',arrayTag,spTag),mArray_sin)
        numpy.savetxt('%s%s_%s_%s.dat'%(mbbDir,mbbFilename+'_cos2',arrayTag,spTag),mArray_cos2)
        numpy.savetxt('%s%s_%s_%s.dat'%(mbbDir,mbbFilename+'_sin2',arrayTag,spTag),mArray_sin2)
        numpy.savetxt('%s%s_%s_%s.dat'%(mbbDir,mbbFilename+'_sincos',arrayTag,spTag),mArray_cossin)
                   
        numpy.savetxt('%sBbl.%s_%s_%s.dat'%(mbbDir,mbbFilename,arrayTag,spTag),Bbl)
        numpy.savetxt('%sBbl.%s_%s_%s.dat'%(mbbDir,mbbFilename+'_cos1',arrayTag,spTag),Bbl_cos)
        numpy.savetxt('%sBbl.%s_%s_%s.dat'%(mbbDir,mbbFilename+'_sin1',arrayTag,spTag),Bbl_sin)
        numpy.savetxt('%sBbl.%s_%s_%s.dat'%(mbbDir,mbbFilename+'_cos2',arrayTag,spTag),Bbl_cos2)
        numpy.savetxt('%sBbl.%s_%s_%s.dat'%(mbbDir,mbbFilename+'_sin2',arrayTag,spTag),Bbl_sin2)
        numpy.savetxt('%sBbl.%s_%s_%s.dat'%(mbbDir,mbbFilename+'_sincos',arrayTag,spTag),Bbl_cossin)


    return mArray


def computeBinMap(iy,ix,powerMaskHolder_11_ly,powerMaskHolder_11_lx,powerMaskShifted_11,powerMaskShifted_12,powerMaskShifted_22,cos_array,sin_array,binMap,binMap0, binMap_cos,binMap_sin,binMap_cos2,binMap_sin2,binMap_cossin,phlx,phly,type):
    
    for i in xrange(len(iy)):
        
        Ly = powerMaskHolder_11_ly[iy[i]]
        Lx = powerMaskHolder_11_lx[ix[i]]
        ang2=numpy.arctan2(Ly,Lx)
        cos=numpy.cos(2*ang2)*cos_array+numpy.sin(2*ang2)*sin_array
        sin=numpy.sin(2*ang2)*cos_array-numpy.cos(2*ang2)*sin_array
        
        trimMap_11= binMap[location[0][i],location[1][i]]*(trimShiftKMap(powerMaskShifted_11,trimAtL,Lx,Ly,phlx,phly,type))
        trimMap_12= binMap[location[0][i],location[1][i]]*(trimShiftKMap(powerMaskShifted_12,trimAtL,Lx,Ly,phlx,phly,type))
        trimMap_22= binMap[location[0][i],location[1][i]]*(trimShiftKMap(powerMaskShifted_22,trimAtL,Lx,Ly,phlx,phly,type))
        
        
        binMap0 += trimMap_11
        binMap_cos += trimMap_12*cos
        binMap_sin += trimMap_12*sin
        binMap_cos2 += trimMap_22*cos**2
        binMap_sin2 += trimMap_22*sin**2
        binMap_cossin += trimMap_22*cos*sin

