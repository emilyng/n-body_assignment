class Potential:
    acceleration = -9.8
    def compute(self, x, y, z, t, mass):
        return (0.0, mass * self.acceleration, 0.0)
