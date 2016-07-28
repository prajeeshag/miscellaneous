#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

dat=np.loadtxt('output.dat')
int=1
print(int)
#plt.plot(dat[:,0], dat[:,1], 'r')
#plt.axis([np.min(dat[:,0]), np.max(dat[:,0]), np.min(dat[:,1]), np.max(dat[:,1])])
#plt.show()

"""
A simple example of an animated plot
"""

fig, ax = plt.subplots()

line, = ax.plot([0, 0], [0, -100], 'ro')

def animate(i):
    line.set_ydata(dat[i,1])  # update the data
    line.set_xdata(0)  # update the data
    return line,

#Init only required for blitting to give a clean slate.
def init():
    line.set_ydata(np.ma.array(x, mask=True))
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(0, dat[:,0].size),
    interval=int, repeat=False, blit=False)

plt.show()
