from pyevtk.hl import imageToVTK
import dicom
import fnmatch
import os
import sys
import numpy as np
#import tifffile 
import glob as gb 
import operator
from scipy.ndimage import zoom

def GenVTK(XCrop, YCrop, XYScale):
    Names=sorted(gb.glob('IM*'))
    #XCrop, YCrop, XYScale =50,100, 25 #%
    Zsize=len(Names)
    Idx = 0
    if Zsize==0:
        Zsize=len(Names) 
    for kat in Names[0:Zsize]:
        DI = dicom.read_file(kat)
        Pix = DI.pixel_array
        Pix=zoom(Pix, XYScale*0.01)
        #Pix2 = Pix
        Xsize, Ysize = Pix.shape
        XPart, YPart = Pix.shape[0]*XCrop/100, Pix.shape[1]*YCrop/100
        #patata
        if kat == Names[0]:
            ###allocate Volume
            Vol=np.zeros([XPart-1, YPart-1, Zsize])
        Vol[:,:,Idx]=Pix[1:XPart,1:YPart]
        Idx+=1
    return Vol;
    #patata
    #imageToVTK("./VTKimage", cellData = {"Scattering" : Vol} )
 #TempFileprefix


#os.chdir(sys.argv[1])
FoldTempFile=sys.argv[1]
TempFileprefix='VTK'
NameMatch='*_OCT'
matches = []
foldNum = []
XCrop = raw_input('enter percentage crop along the X direction:')
YCrop = raw_input('enter percentage crop along the Y direction:')
XYScale = raw_input('enter percentage to scale XY plane:')
XCrop = int(XCrop)
YCrop = int(YCrop)
XYScale = int(XYScale)

for root, dirnames, filenames in os.walk(sys.argv[1]):##to be changed to sys.argv[1]
    for foldname in fnmatch.filter(dirnames, NameMatch):
        pet=os.path.join(root, foldname)
        matches.append(pet)
        petSplt=pet.split('/')
        T=next(i for i in petSplt if i.startswith('PAT'))
        DateAcq=petSplt[petSplt.index(T)+1]
        #DateAcq=pet[pet.index(T)+1]
        val0=next(i for i in petSplt if i.startswith('export'))
        val1=next(i for i in petSplt if i.endswith('_OCT'))
        val1=int(val1.split('_OCT')[0])
        foldNum.append([val0,DateAcq,val1])#int(TMP[len(TMP)-1]))
#patata
list1=sorted(foldNum, key=operator.itemgetter(0, 1, 2))
matchesSort=[]
#patata###This loop converts the elements of the path that were sorted into the path appropriate path
for katL in list1:
    Item=0#for katM in matches:
    while (~np.isnan(Item)):
        katM=matches[Item]
        if (katL[0] in katM and '/'+str(katL[1])+'/' in katM and '/'+str(katL[2])+'_OCT' in katM):
            matchesSort.append(katM)
            Item=np.nan
        Item+=1
#SortPos = np.argsort(foldNum)
TimeP = len(matchesSort)
OoM=str(TimeP).count('')-1
OutFmt="%0"+str(OoM)+"d"
TimeP=0
for katF in matchesSort:
    #CWDst = os.getcwd()
    os.chdir(katF)
    #patata
    VolVTK=GenVTK()
    FoldOut=katF.split('export')[0]
    #os.chdir(CWDst)
    OutNum=OutFmt % (TimeP,)
    imageToVTK(FoldOut+TempFileprefix+"_"+OutNum, cellData = {"Scattering" : VolVTK} )
    TimeP+=1
    print(TempFileprefix+"_"+OutNum)
 
