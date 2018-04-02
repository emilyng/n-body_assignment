from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.animation as animation
import matplotlib.patches as patches
import random

theta = 5
G = 6.67408e-11 #m^3 kg^-1 s^-2
fig1 = plt.figure()
sim = fig1.add_subplot(111, aspect='equal')
fig2 = plt.figure()
ax = plt.axes(projection = '3d')
quadt = fig2.add_subplot(111, aspect='equal')

class Node: 
    children = None
    mass = None
    center_of_mass = None
    bbox = None
    vx = vy = None

def quad_insert(root, x, y, z, m):
    if root.mass is None:   #when the root is empty, add the first particle
        root.mass = m
        root.center_of_mass = np.array([x, y, z])
        return
    elif root.children is None:
        root.children = [None,None,None,None,None,None,None,None]
        old_quadrant = quadrant_of_particle(root, root.center_of_mass[0], 
                                                  root.center_of_mass[1],
                                                  root.center_of_mass[2])
        if root.children[old_quadrant] is None:
            root.children[old_quadrant] = Node()
            root.children[old_quadrant].bbox = quadrant_bbox(root.bbox,old_quadrant)
        quad_insert(root.children[old_quadrant], root.center_of_mass[0], 
                                                 root.center_of_mass[1], 
                                                 root.center_of_mass[2], root.mass)
        new_quadrant = quadrant_of_particle(root, x, y, z)
        if root.children[new_quadrant] is None:
            root.children[new_quadrant] = Node()
            root.children[new_quadrant].bbox = quadrant_bbox(root.bbox,new_quadrant)
        quad_insert(root.children[new_quadrant], x, y, z, m)
        root.center_of_mass[0] = (root.center_of_mass[0]*root.mass + x*m) / (root.mass + m)
        root.center_of_mass[1] = (root.center_of_mass[1]*root.mass + y*m) / (root.mass + m)
        root.center_of_mass[2] = (root.center_of_mass[2]*root.mass + z*m) / (root.mass + m)
        root.mass = root.mass + m
    else:
        new_quadrant = quadrant_of_particle(root, x, y, z)
        if root.children[new_quadrant] is None:
            root.children[new_quadrant] = Node()
            root.children[new_quadrant].bbox = quadrant_bbox(root.bbox, new_quadrant)
        quad_insert(root.children[new_quadrant], x, y, z, m)
        root.center_of_mass[0] = (root.center_of_mass[0]*root.mass + x*m) / (root.mass + m)
        root.center_of_mass[1] = (root.center_of_mass[1]*root.mass + y*m) / (root.mass + m)
        root.center_of_mass[2] = (root.center_of_mass[2]*root.mass + z*m) / (root.mass + m)        
        root.mass = root.mass + m


def display(root):
    if root.mass is None:
        return
    if root.children is not None:
        x = (root.bbox[0] + root.bbox[1]) / 2
        y = (root.bbox[2] + root.bbox[3]) / 2
        width = x-root.bbox[0]
        #plt_node(root.bbox[0], root.bbox[2], width)
        #plt_node(root.bbox[0], y, width)
        #plt_node(x, root.bbox[2], width)
        #plt_node(x, y, width)
        for i in range(4):
            if root.children[i] is not None:
                display(root.children[i])
    else:
        ax.scatter3D(root.center_of_mass[0], root.center_of_mass[1], root.center_of_mass[2])

def integrate(particles):
    bodies = particles
    n = len(bodies)
    t = 0.0
    t_final = 10.0
    dt = 0.1
    while t < t_final:
        particles_force = []
        root = Node()
        root.center_of_mass = []
        root.bbox = find_root_bbox(bodies)
        for i in range(n):
            quad_insert(root, bodies[i].x, bodies[i].y, bodies[i].z, bodies[i].mass)
        for i in range(n):
            total_fx, total_fy, total_fz = compute_force(root, bodies[i].x, bodies[i].y, 
                                                     bodies[i].z, bodies[i].mass)
            particles_force.append((total_fx, total_fy, total_fz))
        for i in range(n):
            fx, fy, fz = particles_force[i]
            bodies[i].Fx += fx / bodies[i].mass * dt
            bodies[i].Fy += fy / bodies[i].mass * dt
            bodies[i].Fz += fz / bodies[i].mass * dt

            bodies[i].x += bodies[i].Fx * dt
            bodies[i].y += bodies[i].Fy * dt
            bodies[i].z += bodies[i].Fz * dt
            sim.scatter(bodies[i].x, bodies[i].y, bodies[i].z)
        t += dt
    display(root)
    ax.scatter3D(root.center_of_mass[0], root.center_of_mass[1], 
                 root.center_of_mass[2], c='red', marker='x')

def compute_force(root, x, y, z, m):
    if root.mass is None:
        return 0, 0, 0
    if root.center_of_mass[0] == x and root.center_of_mass[1] == y and \
       root.center_of_mass[2] == z and root.mass == m:
        return 0, 0, 0
    d = root.bbox[1]-root.bbox[0]
    r = distance(x, y, z, root.center_of_mass[0], 
                          root.center_of_mass[1], 
                          root.center_of_mass[2])
    if d/r < theta or root.children is None:
        return force(m, x, y, z, root.mass, root.center_of_mass[0], 
                                            root.center_of_mass[1],
                                            root.center_of_mass[2])
    else:
        fx = 0.0
        fy = 0.0
        fz = 0.0
        for i in range(8):
            if root.children[i] is not None:
                fx += compute_force(root.children[i], x, y, z, m)[0]
                fy += compute_force(root.children[i], x, y, z, m)[1]
                fz += compute_force(root.children[i], x, y, z, m)[2]
        return fx, fy, fz
    

################################################# SUPPORTING FUNCTION ##############################################################

def force(m, x, y, z, mcm, xcm, ycm, zcm):
    d = distance(x, y, z, xcm, ycm, zcm)
    f = G*m*mcm/(d**2)
    dx = xcm - x
    dy = ycm - y
    dz = zcm - z
    angle = math.atan2(dy, dx, dz)      ##################
    fx = math.cos(angle) * f            ##################
    fy = math.sin(angle) * f            #########force/distance functions need to be fixed!!!!
    fz = math.tan(angle) * f            ##################
    return fx, fy, fz                   ##################

def distance(x1, y1, z1, x2, y2, z2):
    return math.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

'''
def plt_node(x, y, z, width):
    quadt.add_patch(patches.Rectangle((x, y), width, width, fill = False))
'''
'''
def find_root_bbox(particles):
    """ Create a suitable square boundary box for the input particles
    """
    if len(particles) == 0 or len(particles) == 1:
        return None
    x_list = []
    y_list = [] 
    z_list = []
    for particle in particles:
        x_list.append(particle.x)
        y_list.append(particle.y)
        z_list.append(particle.z)
        
    xmin = np.min(x_list)
    xmax = np.max(x_list)
    ymin = np.min(y_list)
    ymax = np.min(y_list)
    zmin = np.min(z_list)
    zmax = np.max(z_list)
        
    for particle in particles:
        if particle.x > xmax:
            xmax = particle.x
        if particle.x < xmin:
            xmin = particle.x
        if particle.y > ymax:
            ymax = particle.y
        if particle.y < ymin:
            ymin = particle.y
        if particle.z > zmax:
            zmax = particle.z
        if particle.z < zmin:
            zmin = particle.z
            
    if xmax - xmin == ymax - ymin == zmax - zmin:
        return xmin, xmax, ymin, ymax, zmin, zmax
    elif xmax - xmin > ymax - ymin:
        return xmin, xmax, ymin, ymax+(xmax-xmin-ymax+ymin)
    else:
        return xmin, xmax+(ymax-ymin-xmax+xmin), ymin, ymax
'''
def find_root_bbox(particles):
    """ Create a suitable square boundary box for the input particles
    """
    if len(particles) == 0 or len(particles) == 1:
        return None
    x_list = []
    y_list = [] 
    z_list = []
    for particle in particles:
        x_list.append(particle.x)
        y_list.append(particle.y)
        z_list.append(particle.z)

    xmin = np.min((np.min(x_list), np.min(y_list), np.min(z_list)))
    ymin = xmin
    zmin = xmin
    xmax = np.max((np.max(x_list), np.max(y_list), np.max(z_list)))
    ymax = xmax
    zmax = xmax
        
    return (xmin, ymin, zmin, xmax, ymax, zmax)

def quadrant_of_particle(root, x, y, z):
    """Return position of quadrant of the particle (x,y,z)
    """
    mid_x = (root.bbox[3] - root.bbox[0])/2
    mid_y = (root.bbox[4] - root.bbox[1])/2
    mid_z = (root.bbox[5] - root.bbox[2])/2
    #midpoint = np.array([mid_x, mid_y, mid_z])   

    if y >= mid_y:
        if x <= mid_x:
            if z >= mid_z:
                return 0
            else:
                return 1
        elif z >= mid_z:
            return 2
        else:
            return 3
    else:
        if x <= mid_x:
            if z >= mid_z:
                return 4
            else:
                return 5
        elif z >= mid_z:
            return 6
        else:
            return 7
    
    root.bbox = np.array(root.bbox)

def quadrant_bbox(bbox, quadrant):
    """Return the coordinate of the quadrant
    """
    x = (root.bbox[3] - root.bbox[0])/2
    y = (root.bbox[4] - root.bbox[1])/2
    z = (root.bbox[5] - root.bbox[2])/2
    
    #Quadrant 0: (xmin, x, y, ymax, z, zmax)
    if quadrant == 0:
        return bbox[0], x, y, bbox[4], z, bbox[5]
    #Quadrant 1: (xmin, x, y, ymax, zmin, z)
    elif quadrant == 1:
        return bbox[0], x, y, bbox[4], bbox[2], z
    #Quadrant 2: (x, xmax, y, ymax, z, zmax)
    elif quadrant == 2:
        return x, bbox[3], y, bbox[4], z, bbox[5]
    #Quadrant 3: (x, xmax, y, ymax, zmin, z)
    elif quadrant == 3:
        return x, bbox[3], y, bbox[4], bbox[2], z
    #Quadrant 4: (xmin, x, ymin, y, z, zmax)
    elif quadrant == 4:
        return bbox[0], x, bbox[1], y, z, bbox[5]
    #Quadrant 5: (xmin, x, ymin, y, zmin, z)
    elif quadrant == 5:
        return bbox[0], x, bbox[1], y, bbox[2], z
    #Quadrant 6: (x, xmax, ymin, y, z, zmax)
    elif quadrant == 6:
        return x, bbox[3], bbox[1], y, z, bbox[5]
    #Quadrant 7: (x, xmax, ymin, y, zmin, z)
    else:
        return x, bbox[3], bbox[1], y, bbox[2], z

'''
def data_from_file(filename, array):
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                continue
            else:
                name,color,m,x,y,vx,vy = line.split(',')
                array.append([name,color,float(m),float(x)*AU,float(y)*AU,float(vx)*1000,float(vy)*1000])
'''

if __name__ == '__main__':
    #filename = ('solar-system.txt')
    #particles = []
    #data_from_file(filename, particles)
    from hundred_particles_list import particles
    root = Node()
    root.center_of_mass = []
    root.bbox = find_root_bbox(particles[:4])
    for particle in particles:
        quad_insert(root, particle.x, particle.y, particle.z, particle.mass)
    print('Boundary box: ',root.bbox)
    print('Total mass: ',root.mass)
    print('Coordinate of center of mass: ',root.center_of_mass)
    plt.scatter(root.center_of_mass[0], root.center_of_mass[1], c='r', marker='x', s=50)
    print('Theta: ', theta)
    integrate(particles)
    plt.show()