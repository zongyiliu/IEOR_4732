# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 11:32:18 2018

@author: Ali
"""

import matplotlib
matplotlib.rcParams['backend'] = "WXAgg"
import matplotlib.pyplot as plt

F = lambda x: np.sin(2*x)

plt.ion()    
x = np.linspace(0, 1, 200)
plt.plot(x, F(x))


for i in range(100):
    if 'ax' in globals(): ax.remove()
    newx = np.random.choice(x, size = 10)
    ax = plt.scatter(newx, F(newx))
    plt.pause(0.05)

plt.ioff()
plt.show()