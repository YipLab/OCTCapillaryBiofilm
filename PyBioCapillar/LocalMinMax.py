import numpy as np

def posMaxArr(ys,w):
    yL=np.roll(ys,w)
    yR=np.roll(ys,-w)
    yT=np.zeros([len(ys),3])
    yT[:,0]=yL
    yT[:,1]=ys
    yT[:,2]=yR
    yTM=yT.max(axis=1)
    args=np.arange(len(ys))
    #len(args)
    OutPut=args[yTM==ys]
    return OutPut
    #yR=np.roll(ys,-w)


def posMinArr(ys,w):
    yL=np.roll(ys,w)
    yR=np.roll(ys,-w)
    yT=np.zeros([len(ys),3])
    yT[:,0]=yL
    yT[:,1]=ys
    yT[:,2]=yR
    yTM=yT.min(axis=1)
    args=np.arange(len(ys))
    #len(args)
    OutPut=args[yTM==ys]
    return OutPut
 
