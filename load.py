import numpy as np
from lib_zjs import draw as vi
from matplotlib import pyplot as plt
step=300
sizeX=101
sizeY=101
#fid = open("tmpBx.dat",'rb')
#Bx=np.fromfile(fid, dtype="float64",count=sizeX*(sizeY-1)*step).reshape((step,sizeY-1,sizeX))
#fid = open("tmpBy.dat",'rb')
#By=np.fromfile(fid, dtype="float64",count=(sizeX-1)*sizeY*step).reshape((step,sizeY,sizeX-1))
#Ez=np.fromfile(fid, dtype="float64",count=(sizeX-1)*(sizeY-1)*step).reshape((sizeY,sizeX,step),order='F')


fid = open("tmpEx.dat",'rb')
Ex=np.fromfile(fid, dtype="float64",count=(sizeX)*(sizeY-1)*step).reshape((step,sizeY-1,sizeX)).transpose()
fid = open("tmpEz.dat",'rb')
Ez=np.fromfile(fid, dtype="float64",count=(sizeX)*(sizeY)*step).reshape((step,sizeY,sizeX)).transpose()

#fid = open("tmpBz.dat",'rb')
#Bz=np.fromfile(fid, dtype="float64",count=(sizeX-1)*(sizeY-1)*step).reshape((sizeY-1,sizeX-1,step),order='F')
