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
        
        self.Fx = 0.0
        self.Fy = 0.0
        self.Fz = 0.0   
        
        self.t = 0.0
        
    '''Basic Euler integration'''   
    def Euler_update(self, Fx, Fy, Fz, dt):
        self.x = self.x + self.vx * dt
        self.y = self.y + self.vy * dt
        self.z = self.z + self.vz * dt
        
        self.vx = self.vx + Fx / self.mass * dt
        self.vy = self.vy + Fy / self.mass * dt
        self.vz = self.vz + Fz / self.mass * dt
        
    '''Semi-implicit Euler integration'''
    def Semi_Euler_update(self, Fx, Fy, Fz, dt):
        self.vx = self.vx + Fx / self.mass * dt
        self.vy = self.vy + Fy / self.mass * dt
        self.vz = self.vz + Fz / self.mass * dt
        
        self.x = self.x + self.vx * dt
        self.y = self.y + self.vy * dt
        self.z = self.z + self.vz * dt
        
    '''Leapfrog integration'''
    def Leapfrog_update(self, Fx, Fy, Fz, dt):
        x_half = self.x + 0.5 * self.vx * dt 
        y_half = self.y + 0.5 * self.vy * dt
        z_half = self.z + 0.5 * self.vz * dt
        
        self.vx = self.vx + x_half * Fx / self.mass * dt
        self.vy = self.vy + y_half * Fy / self.mass * dt
        self.vz = self.vz + z_half * Fz / self.mass * dt
        
        self.x = x_half + 0.5 * self.vx * dt
        self.y = y_half + 0.5 * self.vy * dt
        self.z = z_half + 0.5 * self.vz * dt
    
    
    '''[force] METHOD used on particles / Direct calculation of forces'''
    '''Direct calculation of forces'''
    def pairwise_force(self, particle):
        G = 6.67e-4
        
        r2 = (self.x - particle.x)**2.0 + \
             (self.y - particle.y)**2.0 + \
             (self.z - particle.z)**2.0
        F_mag = -(G * particle.mass * self.mass)/r2
        F_x = (self.x - particle.x)/r2**0.5 * F_mag
        F_y = (self.y - particle.y)/r2**0.5 * F_mag
        F_z = (self.z - particle.z)/r2**0.5 * F_mag
        return (F_x, F_y, F_z)
