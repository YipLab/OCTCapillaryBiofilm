import numpy as np

def posMinArr(ys,w):
    yL=np.roll(ys,w)
    yR=np.roll(ys,-w)
    yT=zeros(len(ys),3)
    yT[:,0]=yL
    yT[:,1]=ys
    yT[:,2]=yR
    yT.max(axis=1)
    args=np.arange(len(ys))
    OutPut=args[yT==ys]
    return OutPut
    #yR=np.roll(ys,-w)
