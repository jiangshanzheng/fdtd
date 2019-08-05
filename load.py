import numpy as np
mod=1;
step=300
sizeX=128
sizeY=128
#fid = open("tmpBx.dat",'rb')
#Bx=np.fromfile(fid, dtype="float64",count=sizeX*(sizeY-1)*step).reshape((step,sizeY-1,sizeX))
#fid = open("tmpBy.dat",'rb')
#By=np.fromfile(fid, dtype="float64",count=(sizeX-1)*sizeY*step).reshape((step,sizeY,sizeX-1))

if 0==mod:
    fid = open("tmpEx.dat",'rb')
    Ex=np.fromfile(fid, dtype="float64",count=(sizeX-1)*(sizeY)*step).reshape((step,sizeX-1,sizeY)).transpose()
    fid = open("tmpEy.dat",'rb')
    Ey=np.fromfile(fid, dtype="float64",count=(sizeX)*(sizeY-1)*step).reshape((step,sizeX,sizeY-1)).transpose(1,2,0)
#Ey=np.fromfile(fid, dtype="float64",count=(sizeX)*(sizeY-1)*step).reshape((sizeX,sizeY-1,step),order='F')
else:
    fid = open("tmpEz.dat",'rb')
    Ez=np.fromfile(fid, dtype="float64",count=(sizeX)*(sizeY)*step).reshape((step,sizeY,sizeX)).transpose()

#fid = open("tmpBz.dat",'rb')
#Bz=np.fromfile(fid, dtype="float64",count=(sizeX-1)*(sizeY-1)*step).reshape((sizeY-1,sizeX-1,step),order='F')
