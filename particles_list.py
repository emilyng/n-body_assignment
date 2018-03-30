#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 19:04:29 2018

@author: emilyng
"""

import h5py

from particle import Particle

'''Read into hdf5 file with initial conidtions'''
'''file names 'gaussian.h5' located in working directory'''

hf = h5py.File("gaussian.h5", 'r')

'''retrieve particles positions, velocities, and masses from hdf5 file'''
positions = hf['particle_positions'][:]
velocities = hf['particle_velocities'][:]
masses = hf['particle_masses'][:]

'''Create list of particle objects from initial conditions'''
particles = []
for i in range(len(positions)):
    particles.append(Particle(positions[i][0], positions[i][1], positions[i][2],
                              velocities[i][0], velocities[i][1], velocities[i][2],
                              masses[i]))


 