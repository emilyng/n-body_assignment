class System:
    """ Creates a '''system''' of particles with a background potential"""
    def __init__(self, potential, particles):
        self.potential = potential
        self.particles = particles
        self.t = 0.0
        
    def force_calculation(self, body_index):
        N = len(self.particles) 
        force = np.zeros(3)
        forces = []
        
        target_body = self.particles[body_index]
        for i, external_body in enumerate(self.particles):
            if i != body_index:
                    force += self.particles[i].pairwise_force(self.particles[body_index])
            forces.append(list(force))
            
        return forces
    
    def force_sum(self):
        N = len(self.particles) 
        
        for body_index in range(N):
            forces = self.force_calculation(body_index)
    
        for i in range(N):
            self.particles[i].Fx = forces[i][0]
            self.particles[i].Fy = forces[i][1]
            self.particles[i].Fz = forces[i][2]
            
        
    def Euler_update(self, dt):
        for particle in self.particles:
            Fx_b, Fy_b, Fz_b = self.potential.compute(
                particle.x, particle.y, particle.z, self.t, particle.mass)
        
        for i, particle in enumerate(self.particles):
            self.force_sum()
            
            self.particles[i].Fx = particle.Fx + Fx_b
            self.particles[i].Fy = particle.Fy + Fy_b
            self.particles[i].Fz = particle.Fz + Fz_b
            
            Fx = particle.Fx
            Fy = particle.Fy
            Fz = particle.Fz
            particle.Euler_update(Fx, Fy, Fz, dt)
        self.t += dt
        
    def Semi_Euler_update(self, dt):
        for particle in self.particles:
            Fx_b, Fy_b, Fz_b = self.potential.compute(
                particle.x, particle.y, particle.z, self.t, particle.mass)
        
        for i, particle in enumerate(self.particles):
            particle.Fx = particle.Fx + Fx_b
            particle.Fy = particle.Fy + Fy_b
            particle.Fz = particle.Fz + Fz_b
            
            Fx = particle.Fx
            Fy = particle.Fy
            Fz = particle.Fz
            particle.Semi_Euler_update(Fx, Fy, Fz, dt)
        self.t += dt
        
    def Leapfrog_update(self, dt):
        for particle in self.particles:
            Fx_b, Fy_b, Fz_b = self.potential.compute(
                particle.x, particle.y, particle.z, self.t, particle.mass)
        
        for particle in self.particles:
            particle.Fx = particle.Fx + Fx_b
            particle.Fy = particle.Fy + Fy_b
            particle.Fz = particle.Fz + Fz_b
            
            Fx = particle.Fx
            Fy = particle.Fy
            Fz = particle.Fz
            particle.Leapfrog_update(Fx, Fy, Fz, dt)
        self.t += dt         
