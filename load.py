import numpy as np
step=2500
size=200
fid = open("tmpE",'rb')
E=np.fromfile(fid, dtype="float64",count=size*step).reshape((step,size))
fid = open("tmpB",'rb')
B=np.fromfile(fid, dtype="float64",count=size*step).reshape((step,size))
