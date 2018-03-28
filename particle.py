class Particle:
    """
    Creates a '''Particle''' object with x-coord, y-coord, z-coord, x-vel, 
        y-vel,z-vel, and mass attached to it
    """
    def __init__(self, x, y, z, vx, vy, vz, mass):
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.mass = mass
