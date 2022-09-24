import sys
import numpy as np
from scipy import signal
energy_row=sys.argv[1]
num=np.fromstring(energy_row,dtype=float,sep=',')
local_max_num=signal.argrelextrema(num,np.greater)
local_max=' '.join(str(i) for i in local_max_num[0])
print(local_max)