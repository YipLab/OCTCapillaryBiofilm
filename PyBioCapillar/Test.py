import EllipseFitter as ef
import dicom
import glob as gb 
import numpy as np
import scipy.signal as ss 
from skimage import filters, measure, segmentation
from scipy.ndimage import find_objects
import os
import matplotlib.pyplot as pp
from skimage.transform import hough_ellipse
import skimage.morphology as morph 
import skimage.measure as meas
import scipy
from LocalMinMax import posMinArr,posMaxArr
from Smooth import smooth
from Centerization import Centress,MinMaxProj
import LogPolarTransform
import tifffile

##from LogPolarTransform import logpolar_fancy
 
###FROM: http://www.lfd.uci.edu/~gohlke/code/imreg.py.html
##from __future__ import division, print_function


##import numpy
##from numpy.fft import fft2, ifft2, fftshi
#from skimage.draw import ellipse_perimeter
#from skimage.filters import roberts, sobel, scharr, prewitt
#from scipy.ndimage.filters import gaussian_filter
#import javabridge as jv
#import bioformats as bf
#jv.start_vm(class_path=bf.JARS)

##https://pyscience.wordpress.com/2014/09/08/dicom-in-python-importing-medical-image-data-into-numpy-with-pydicom-and-vtk/
##http://scikit-image.org/docs/dev/auto_examples/plot_circular_elliptical_hough_transform.html
##http://docs.opencv.org/master/dd/d49/tutorial_py_contour_features.html#gsc.tab=0

FileOnFold=sorted(gb.glob('IM*'), key=os.path.getmtime)##gb.glob('IM*',)
#import matplotlib.pyplot as pp
#import cv2
AvgWin=5
ExtSrch=500
MaxHough=20
HistS=100
Padd=100
FileNum=len(FileOnFold)
RefDs = dicom.read_file(FileOnFold[0])
Img=RefDs.pixel_array
Stack = np.zeros([Img.shape[0],Img.shape[1],FileNum])
Count=0
for kat in FileOnFold:
    RefDs = dicom.read_file(kat)##rdr = bf.ImageReader(kat, perform_init=True)
    Stack[:,:,Count]=RefDs.pixel_array
    Count+=1
    #ConstPixelDims = (int(RefDs.Rows), int(RefDs.Columns))
    # Load spacing values (in mm)
    #ConstPixelSpacing = (float(RefDs.PixelSpacing[0]), float(RefDs.PixelSpacing[1]), float(RefDs.SliceThickness))

ZProj=Stack.mean(axis=2)
#[ZProjLog,Base]=logpolar(ZProj)


ZPLog2=LogPolarTransform.logpolar_fancy(ZProj,ZProj.shape[0]/2,ZProj.shape[1]/2)
MinPos,MaxPos=MinMaxProj(ZPLog2)
RefDiff0=MaxPos-MinPos
Rh,Rv,SegmntSMPadd1=Centress(ZPLog2)
ZProjT=ZProj
ZP1=np.roll(ZProjT,Rv,axis=0)
ZP2=np.roll(ZP1,Rh,axis=1)
ZP2Log=LogPolarTransform.logpolar_fancy(ZP2,ZP2.shape[0]/2,ZP2.shape[1]/2)
MinPos1,MaxPos1=MinMaxProj(ZP2Log)
RefDiff=MaxPos1-MinPos1
Rh=Rh/2;Rv=Rv/2
#if RefDiff<RefDiff0:
Choose=True
while (np.abs(Rh)>1 and np.abs(Rv)>1):
    if Choose:
        RefDiff0=RefDiff
        ZP1=np.roll(ZP2,Rv,axis=0)
        ZP2=np.roll(ZP1,Rh,axis=1)
        ZP2Log=LogPolarTransform.logpolar_fancy(ZP2,ZP2.shape[0]/2,ZP2.shape[1]/2)
        MinPos1,MaxPos1=MinMaxProj(ZP2Log)
        RefDiff=MaxPos1-MinPos1
    else:
        RefDiff0=RefDiff
        ZP1=np.roll(ZP2,-Rv,axis=0)
        ZP2=np.roll(ZP1,-Rh,axis=1)
        ZP2Log=LogPolarTransform.logpolar_fancy(ZP2,ZP2.shape[0]/2,ZP2.shape[1]/2)
        MinPos1,MaxPos1=MinMaxProj(ZP2Log)
        RefDiff=MaxPos1-MinPos1
    ##RefDiff0-RefDiff
    if ((-1)**Choose*(RefDiff0-RefDiff))>=0:
        Rh=Rh/2;Rv=Rv/2
        Choose=False
        ##Rh,Rv
        




Rh,Rv,SegmntSMPadd2=Centress(ZP2Log)
ZProjT=ZP2
ZP1=np.roll(ZProjT,Rv,axis=0)
ZP2=np.roll(ZP1,Rh,axis=1)
ZP2Log=LogPolarTransform.logpolar_fancy(ZP2,ZP2.shape[0]/2,ZP2.shape[1]/2)
Rh,Rv,SegmntSMPadd3=Centress(ZP2Log)
ZProjT=ZP2
ZP1=np.roll(ZProjT,Rv,axis=0)
ZP2=np.roll(ZP1,Rh,axis=1)
ZP2Log=LogPolarTransform.logpolar_fancy(ZP2,ZP2.shape[0]/2,ZP2.shape[1]/2)
Rh,Rv,SegmntSMPadd4=Centress(ZP2Log)
ZProjT=ZP2
ZP1=np.roll(ZProjT,Rv,axis=0)
ZP2=np.roll(ZP1,Rh,axis=1)
ZP2Log=LogPolarTransform.logpolar_fancy(ZP2,ZP2.shape[0]/2,ZP2.shape[1]/2)
Rh,Rv,SegmntSMPadd5=Centress(ZP2Log)
ZProjT=ZP2
ZP1=np.roll(ZProjT,Rv,axis=0)
ZP2=np.roll(ZP1,Rh,axis=1)
ZP2Log=LogPolarTransform.logpolar_fancy(ZP2,ZP2.shape[0]/2,ZP2.shape[1]/2)
StackLogPol=np.zeros(Stack.shape)
for kat in np.arange(Stack.shape[2]):
    Slice=Stack[:,:,kat]
    ZP1=np.roll(Slice,Rv,axis=0)
    ZP2=np.roll(ZP1,Rh,axis=1)
    ZP2Log=LogPolarTransform.logpolar_fancy(ZP2,ZP2.shape[0]/2,ZP2.shape[1]/2)
    StackLogPol[0:ZP2Log.shape[0]-ZP2Log.shape[0]*3/4,:,kat]=ZP2Log[ZP2Log.shape[0]*3/4:ZP2Log.shape[0],:]
    
STLP=StackLogPol.mean(axis=1)
STLP1=STLP.max(axis=1)
STLP1Pos=np.where(STLP1>STLP1.mean())
MaxPos=STLP1Pos[0].max()
MaxPos=MaxPos+MaxPos*0.25
DataSave=StackLogPol[0:MaxPos,:,:]
DataSave=DataSave.astype(np.int32)

#pp.matshow(ZP2Log)
#pp.show()
FileOut='Straight_'+FileOnFold[0][0:2]
for kat in np.arange(DataSave.shape[1]):
    OutPutSingleFile = FileOut + str(kat) + '.tif'
    tifffile.imsave(OutPutSingleFile,DataSave[:,kat,:])


with tifffile.TiffWriter(OutPutSingleFile, bigtiff=True) as tif:
    for i in range(DataSave.shape[0]):
        tif.save(DataSave[i], compress=0)

with tifffile.TiffWriter(OutPutSingleFile, bigtiff=True) as tif:
    tif.imsave(DataSave, compress=0)

patata
ZPLog2Mn=ZProjLog2.mean(axis=1)
Mx=ZPLog2Mn.max()
Md=np.median(ZPLog2Mn)
RPos=np.arange(len(ZPLog2Mn))
RPosNaN=RPos+np.nan
ValidRPos=np.where(ZPLog2Mn>((Mx+Md)/2),RPos,RPosNaN)
ValidRPos=ValidRPos[~np.isnan(ValidRPos)].astype(int)
MinMax=ValidRPos.min()-len(ZPLog2Mn)/20
MaxMax=min(ValidRPos.max()+len(ZPLog2Mn)/20,len(ZPLog2Mn))
##pp.plot(RPos[MinMax:MaxMax],ZPLog2Mn[MinMax:MaxMax],'xb')


Segmnt=ZProjLog2[MinMax:MaxMax,:]
SegmntMx=Segmnt.argmax(axis=0)
SegmntSM=smooth(Segmnt.argmax(axis=0),window_len=11)#Segmnt.argmax(axis=0)
MaxPadd=SegmntSM[0:2*ExtSrch]
MinPadd=SegmntSM[len(SegmntSM)-2*ExtSrch:len(SegmntSM)]
SegmntSMPadd=np.zeros(len(SegmntSM)+4*ExtSrch)
PosMaxCtr=SegmntSM.argmax()+2*ExtSrch
PosMinCtr=SegmntSM.argmin()+2*ExtSrch
SegmntSMPadd[0:2*ExtSrch]=MinPadd
SegmntSMPadd[2*ExtSrch:len(SegmntSM)+2*ExtSrch]=SegmntSM
SegmntSMPadd[len(SegmntSM)+2*ExtSrch:len(SegmntSMPadd)]=MaxPadd
#pp.plot(SegmntSMPadd)
MinPos=posMinArr(SegmntSMPadd,750)
MinCut=np.where(np.diff(MinPos)!=1)[0][0]
MinPos1=MinPos[0:MinCut+1].mean()
MinPos2=MinPos[MinCut+1:len(MinPos)].mean()
if MinPos1>(len(SegmntSMPadd)-MinPos2):
    MinPos=MinPos1
else:
    MinPos=MinPos2
    
MinPosReal=MinPos-2*ExtSrch

MaxPos=posMaxArr(SegmntSMPadd,750)
MaxCut=np.where(np.diff(MaxPos)!=1)[0][0]
MaxPos1=MaxPos[0:MaxCut+1].mean()
MaxPos2=MaxPos[MaxCut+1:len(MaxPos)].mean()
if MaxPos1>(len(SegmntSMPadd)-MaxPos2):
    MaxPos=MaxPos1
else:
    MaxPos=MaxPos2

MaxPosReal=MaxPos-2*ExtSrch
MaxAng=MaxPosReal*2*np.pi/len(SegmntSM)
MinAng=MinPosReal*2*np.pi/len(SegmntSM)
AvgAng=(MaxAng-np.pi+MinAng)/2##180=np.pi convert both angles to the same reference
RadShift=(SegmntMx[MaxPosReal]-SegmntMx[MinPosReal])/2
#SegmntMx[MinPosReal]
#SegmntMx[MaxPosReal]
#pp.plot(MaxPos,SegmntSMPadd[MaxPos],'+b')
if AvgAng<(np.pi/2):
    Rh=RadShift*np.cos(AvgAng)
    Rv=RadShift*np.sin(AvgAng)
elif AvgAng<np.pi:
    Rh=-RadShift*np.cos(np.pi-AvgAng)
    Rv=RadShift*np.sin(np.pi-AvgAng)
elif AvgAng<(3*np.pi/4):
    Rh=-RadShift*np.cos(0.75*np.pi-AvgAng)
    Rv=-RadShift*np.sin(0.75*np.pi-AvgAng)
else:
    Rh=RadShift*np.cos(2*np.piAvgAng)
    Rv=-RadShift*np.sin(2*np.piAvgAng)
Rh=int(np.round(Rh))
Rv=int(np.round(Rv))
ZP1=np.roll(ZProj,Rv,axis=0)
ZP2=np.roll(ZP1,Rh,axis=1)
ZP2Log=LogPolarTransform.logpolar_fancy(ZP2,ZP2.shape[0]/2,ZP2.shape[1]/2)



ptata






FitXs=np.arange(PosMinCtr-ExtSrch,PosMinCtr+ExtSrch)
FitPnts=SegmntSMPadd[FitXs]

MinPos=posMinArr(FitPnts,1)
FitPnts[MinPos].argmin()

FitXs=np.arange(PosMaxCtr-ExtSrch,PosMaxCtr+ExtSrch)
FitPnts=SegmntSMPadd[FitXs]

MaxPos=posMaxArr(FitPnts,1)
FitPnts[MaxPos].argmax()



pp.plot(FitPnts,'ro')
pp.plot(OutPut,FitPnts[OutPut])
f=interp1d(FitXs,FitPnts, kind='cubic')

NewXs=np.arange(PosMinCtr-ExtSrch,PosMinCtr+ExtSrch,10)
#[a,b,c]=np.polyfit(FitXs,FitPnts,2)
#Xmax=

FitPnts=SegmntSMPadd[PosMaxCtr-ExtSrch:PosMaxCtr+ExtSrch]



Col=ss.medfilt(ZProjLog.max(axis=0),25)
MaxPos=Col.argmax()
scipy.misc.imsave('ZProjLog.jpg', ZProjLog)
scipy.misc.imsave('ZProj.jpg', ZProj)
#HoughRadii=np.arange(MaxPos-MaxHough,MaxPos+MaxHough,5)
#hough_res = hough_circle(ZProj, HoughRadii)

M=meas.moments(ZProj.astype('uint8'))
cx = M[0, 1] / M[0, 0]
cy = M[1, 0] / M[0, 0]
[cxc,cyc]=ZProj.shape
yCorr=int(cy-cyc/2)
xCorr=int(cx-cxc/2)
ZProj1=np.roll(ZProj,-xCorr,axis=0)
ZProj2=np.roll(ZProj1,-yCorr,axis=1)
pp.matshow(ZProj2)




Col=ss.medfilt(Img.max(axis=0),5)
Row=ss.medfilt(Img.max(axis=1),5)
MR=np.median(Row)
MC=np.median(Col)
TFC=Col>(MC-MC/10)
TFR=Row>(MR-MR/10)
CP=np.where(TFC==True)
RP=np.where(TFR==True)
TopCut=np.max([np.min(RP)-Padd,0])
BotCut=np.min([np.max(RP)+Padd,len(Row)])
LeftCut=np.max([np.min(CP)-Padd,0])
RightCut=np.min([np.max(CP)+Padd,len(Col)])
StackC=Stack[TopCut:BotCut,LeftCut:RightCut,:]
#ZProj=StackC.max(axis=2)
RefDs = dicom.read_file(FileOnFold[FileNum-1])
Img=RefDs.pixel_array
Img=Img[TopCut:BotCut,LeftCut:RightCut]
ZProjT=filters.threshold_adaptive(filters.gaussian_filter(Img,25),1500)
all_labels = measure.label(ZProjT)
slices=find_objects(all_labels)
NumObj=len(slices)
Areas=np.zeros(NumObj)
for kat in range(NumObj):
    loc=all_labels==kat   
    prop = measure.regionprops(loc)  
    val=[prope.area for prope in prop]##[prope.equivalent_diameter for prope in prop]
    Areas[kat]=val[0] 
    ##pp.matshow(loc)
    ##pp.show()
    
Area=np.fliplr([np.sort(Areas)])[0][2]
Ring=np.where(Areas==Area)
loc=all_labels*0
loc[all_labels==Ring]=255
M=meas.moments(loc.astype('uint8'))

#pp.matshow(loc)
#edges = canny(ImgC, sigma=2.0, low_threshold=0.55, high_threshold=0.8)
    

#M=cv2.moments(ImgC)
#cx = int(M['m10']/M['m00'])
#cy = int(M['m01']/M['m00'])
#Middle=(BotCut-TopCut)/2
    
       #[Count, Bins]=np.histogram(Col,HistS)
    #BinW=np.diff(Bins)/2
    #B_c=Bins[np.arange(HistS)]+BinW
    
    #edges = canny(Img, sigma=3, low_threshold=10, high_threshold=50)
    ##result = hough_ellipse(edges, accuracy=20, threshold=250,min_size=100, max_size=120)
