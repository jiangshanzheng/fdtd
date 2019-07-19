import numpy as np
from lib_zjs import draw as vi
from matplotlib import pyplot as plt
step=2500
size=200
fid = open("tmpE",'rb')
E=np.fromfile(fid, dtype="float64",count=size*step).reshape((step,size))
fid = open("tmpB",'rb')
B=np.fromfile(fid, dtype="float64",count=size*step).reshape((step,size))
A=plt.subplots();plt.ylim(-1,1);vi.cartoon(A,E[:,:].transpose(),0.5,'1d',interval=100)