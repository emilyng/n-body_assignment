#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 19:46:51 2018

@author: emilyng
"""

import numpy as np

from particles_list import particles

            
'''direct calculation (summation of forces)'''
N = len(particles) 

def force_sum(particles):
    N = len(particles) 
    force = np.zeros(3)
    forces = []

    for i in range(N):
        for j in range(N):
            if i != j:
                force += particles[i].pairwise_force(particles[j])
        forces.append(list(force))
    
    return forces

forces = force_sum(particles)

for i in range(N):
    particles[i].Fx = forces[i][0]
    particles[i].Fy = forces[i][1]
    particles[i].Fz = forces[i][2]

################################################################################
################################################################################
################################################################################
    
#class Euler_integrate(particles, )
    
from itertools import repeat
xs = [[] for i in repeat(None, N)]
ys = [[] for i in repeat(None, N)]
zs = [[] for i in repeat(None, N)]

ts = [[] for i in repeat(None, N)]

t_final = 1.0
dt = 0.1

'''Basic Euler integration / Constant potential'''
while t_final > particles[0].t:
    for i, particle in enumerate(particles):
        xs[i].append(particle.x)
        ys[i].append(particle.y)
        zs[i].append(particle.z)
        ts[i].append(particles[0].t)
    
    particles[0].update(dt)
    
    forces = force_sum(particles)

    for i in range(N):
        particles[i].Fx = forces[i][0]
        particles[i].Fy = forces[i][1]
        particles[i].Fz = forces[i][2]
        
'''Plotting the movement of particles'''
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(xs[0], ys[0], zs[0], '.-', label = 'particle 1')
ax.plot(xs[1], ys[1], zs[1], '.-', label = 'particle 2')
ax.plot(xs[2], ys[2], zs[2], '.-', label = 'particle 3')
ax.plot(xs[3], ys[3], zs[3], '.-', label = 'particle 4')
ax.plot(xs[4], ys[4], zs[4], '.-', label = 'particle 5')

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Basic Euler Integration Movement")

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
