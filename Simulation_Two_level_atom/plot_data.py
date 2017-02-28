from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors, ticker, cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

time = []
f1_real = []
f1_imag = []

f2_real = []
f2_imag = []




for line in open ("data_abs.csv"):
    time.append(float(line.split(',')[0]))
    f1_real.append(float(line.split(',')[1]))
    f1_imag.append(float(line.split(',')[2]))
    f2_real.append(float(line.split(',')[3]))
    f2_imag.append(float(line.split(',')[4]))
    
# grid the plot area

gs = gridspec.GridSpec(2, 2)

print(len(time), len(f1_real))


ax1 = plt.subplot(gs[0,0])
plt.plot(time,f1_real)
#plt.xlim(fre_min, fre_max)
ax1.set_ylabel('Real Part')
ax1.set_xlabel('$y\' = x * y$')


ax2 = plt.subplot(gs[0,1])
plt.plot(time,f1_imag)
#plt.xlim(fre_min, fre_max)
ax2.set_ylabel('Imag Part')
ax2.set_xlabel('$y\' = x * y$')

ax3 = plt.subplot(gs[1,0])
plt.plot(time,f2_real)
#plt.xlim(fre_min, fre_max)
ax3.set_ylabel('Real Part')
ax3.set_xlabel('$y\' =  y$')


ax4 = plt.subplot(gs[1,1])
plt.plot(time,f2_imag)
#plt.xlim(fre_min, fre_max)
ax4.set_ylabel('Imag Part')
ax4.set_xlabel('$y\' = y$')



plt.show()
