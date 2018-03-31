class Potential:
    def __init__(self, acceleration):
        self.acceleration = acceleration
    def compute(self, x, y, z, t, mass):
        return (mass * self.acceleration, mass * self.acceleration, mass * self.acceleration)
