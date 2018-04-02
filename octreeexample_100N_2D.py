from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import random

theta = 5
AU = (149.6e6 * 1000)     # 149.6 million km, in meters.
G = 6.67408e-11 #m^3 kg^-1 s^-2
fig1 = plt.figure()
sim = fig1.add_subplot(111, aspect='equal')
fig2 = plt.figure()
quadt = fig2.add_subplot(111, aspect='equal')

class Node: 
    children = None
    mass = None
    center_of_mass = None
    bbox = None
    vx = vy = None

def quad_insert(root, x, y, m):
    if root.mass is None:   #when the root is empty, add the first particle
        root.mass = m
        root.center_of_mass = [x,y]
        return
    elif root.children is None:
        root.children = [None,None,None,None]
        old_quadrant = quadrant_of_particle(root.bbox, root.center_of_mass[0], root.center_of_mass[1])
        if root.children[old_quadrant] is None:
            root.children[old_quadrant] = Node()
            root.children[old_quadrant].bbox = quadrant_bbox(root.bbox,old_quadrant)
        quad_insert(root.children[old_quadrant], root.center_of_mass[0], root.center_of_mass[1], root.mass)
        new_quadrant = quadrant_of_particle(root.bbox, x, y)
        if root.children[new_quadrant] is None:
            root.children[new_quadrant] = Node()
            root.children[new_quadrant].bbox = quadrant_bbox(root.bbox,new_quadrant)
        quad_insert(root.children[new_quadrant], x, y, m)
        root.center_of_mass[0] = (root.center_of_mass[0]*root.mass + x*m) / (root.mass + m)
        root.center_of_mass[1] = (root.center_of_mass[1]*root.mass + y*m) / (root.mass + m)
        root.mass = root.mass + m
    else:
        new_quadrant = quadrant_of_particle(root.bbox, x, y)
        if root.children[new_quadrant] is None:
            root.children[new_quadrant] = Node()
            root.children[new_quadrant].bbox = quadrant_bbox(root.bbox, new_quadrant)
        quad_insert(root.children[new_quadrant], x, y, m)
        root.center_of_mass[0] = (root.center_of_mass[0]*root.mass + x*m) / (root.mass + m)
        root.center_of_mass[1] = (root.center_of_mass[1]*root.mass + y*m) / (root.mass + m)
        root.mass = root.mass + m

def display(root):
    if root.mass is None:
        return
    if root.children is not None:
        x = (root.bbox[0] + root.bbox[1]) / 2
        y = (root.bbox[2] + root.bbox[3]) / 2
        width = x-root.bbox[0]
        plt_node(root.bbox[0], root.bbox[2], width)
        plt_node(root.bbox[0], y, width)
        plt_node(x, root.bbox[2], width)
        plt_node(x, y, width)
        for i in range(4):
            if root.children[i] is not None:
                display(root.children[i])
    else:
        quadt.scatter(root.center_of_mass[0], root.center_of_mass[1])

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
            quad_insert(root, bodies[i].x, bodies[i].y, bodies[i].mass)
        for i in range(n):
            total_fx, total_fy = compute_force(root,bodies[i].x,bodies[i].y,bodies[i].mass)
            particles_force.append((total_fx, total_fy))
        for i in range(n):
            fx, fy = particles_force[i]
            bodies[i].Fx += fx / bodies[i].mass * dt
            bodies[i].Fy += fy / bodies[i].mass * dt

            bodies[i].x += bodies[i].Fx * dt
            bodies[i].y += bodies[i].Fy * dt
            sim.scatter(bodies[i].x, bodies[i].y)
        t += dt
    display(root)
    quadt.scatter(root.center_of_mass[0], root.center_of_mass[1], c='red', marker='x')

def compute_force(root,x,y,m):
    if root.mass is None:
        return 0, 0
    if root.center_of_mass[0] == x and root.center_of_mass[1] == y and root.mass == m:
        return 0, 0
    d = root.bbox[1]-root.bbox[0]
    r = distance(x,y, root.center_of_mass[0], root.center_of_mass[1])
    if d/r < theta or root.children is None:
        return force(m, x, y, root.mass, root.center_of_mass[0], root.center_of_mass[1])
    else:
        fx = 0.0
        fy = 0.0
        for i in range(4):
            if root.children[i] is not None:
                fx += compute_force(root.children[i],x,y,m)[0]
                fy += compute_force(root.children[i],x,y,m)[1]
        return fx, fy

################################################# SUPPORTING FUNCTION ##############################################################

def force(m, x, y, mcm, xcm, ycm):
    d = distance(x, y, xcm, ycm)
    f = G*m*mcm/(d**2)
    dx = xcm - x
    dy = ycm - y
    angle = math.atan2(dy, dx)
    fx = math.cos(angle) * f
    fy = math.sin(angle) * f
    return fx, fy

def distance(x1, y1, x2, y2):
    return math.sqrt((x2-x1)**2+(y2-y1)**2)

def plt_node(x, y, width):
    quadt.add_patch(patches.Rectangle((x, y), width, width, fill = False))

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
    

def quadrant_of_particle(bbox, x, y):
    """Return position of quadrant of the particle (x,y)
    """
    if y >= (bbox[3] + bbox[2])/2:
        if x <= (bbox[1] + bbox[0])/2:
            return 0
        else:
            return 1
    else:
        if x >= (bbox[1] + bbox[0])/2:
            return 2
        else:
            return 3

def quadrant_bbox(bbox,quadrant):
    """Return the coordinate of the quadrant
    """
    x = (bbox[0] + bbox[1]) / 2
    y = (bbox[2] + bbox[3]) / 2
    #Quadrant 0: (xmin, x, y, ymax)
    if quadrant == 0:
        return bbox[0], x, y, bbox[3]
    #Quadrant 1: (x, xmax, y, ymax)
    elif quadrant == 1:
        return x, bbox[1], y, bbox[3]
    #Quadrant 2: (x, xmax, ymin, y)
    elif quadrant == 2:
        return x, bbox[1], bbox[2], y
    #Quadrant 3: (xmin, x, ymin, y)
    elif quadrant == 3:
        return bbox[0], x, bbox[2], y

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
    root.bbox = find_root_bbox(particles)
    for particle in particles:
        quad_insert(root, particle.x, particle.y, particle.mass)
    print('Boundary box: ',root.bbox)
    print('Total mass: ',root.mass)
    print('Coordinate of center of mass: ',root.center_of_mass)
    plt.scatter(root.center_of_mass[0], root.center_of_mass[1], c='r', marker='x', s=50)
    print('Theta: ', theta)
    integrate(particles)
    plt.show()