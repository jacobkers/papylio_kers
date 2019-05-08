# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:05:44 2019

@author: mwdocter
"""

import matplotlib.pyplot as plt
import numpy
a = numpy.arange(int(1e3)) # for testing can use 1e7
# large so you can easily see the memory footprint on the system monitor.
fig = plt.figure()
ax  = fig.add_subplot(1, 2, 1)
pnt1 = ax.plot(a,'b-*') # this uses up an additional 230 Mb of memory.
fig.canvas.draw()
fig.canvas.flush_events()
# can I get the memory back?
l = pnt1[0]
l.remove()
del l
del lines
# not releasing memory
ax.cla() 