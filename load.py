import numpy as np
from lib_zjs import draw as vi
step=300
sizeX=32
sizeY=512
#fid = open("tmpBx.dat",'rb')
#Bx=np.fromfile(fid, dtype="float64",count=sizeX*(sizeY-1)*step).reshape((step,sizeY-1,sizeX))
#fid = open("tmpBy.dat",'rb')
#By=np.fromfile(fid, dtype="float64",count=(sizeX-1)*sizeY*step).reshape((step,sizeY,sizeX-1))


fid = open("tmpEx.dat",'rb')
Ex=np.fromfile(fid, dtype="float64",count=(sizeX-1)*(sizeY)*step).reshape((step,sizeX-1,sizeY)).transpose()


fid = open("tmpEy.dat",'rb')
Ey=np.fromfile(fid, dtype="float64",count=(sizeX)*(sizeY-1)*step).reshape((step,sizeX,sizeY-1)).transpose(1,2,0)
#Ey=np.fromfile(fid, dtype="float64",count=(sizeX)*(sizeY-1)*step).reshape((sizeX,sizeY-1,step),order='F')

#fid = open("tmpEz.dat",'rb')
#Ez=np.fromfile(fid, dtype="float64",count=(sizeX)*(sizeY)*step).reshape((step,sizeY,sizeX)).transpose()

#fid = open("tmpBz.dat",'rb')
#Bz=np.fromfile(fid, dtype="float64",count=(sizeX-1)*(sizeY-1)*step).reshape((sizeY-1,sizeX-1,step),order='F')
