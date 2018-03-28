
from doit.tools import run_once
import h5py
import numpy as np
import matplotlib.pyplot as plt

def task_generate_gaussian():
    N = 32**3
    seed = 0x4d3d3d3
    fn = "gaussian.h5"
    def _generate():
        np.random.seed(seed)
        pos = np.random.normal(loc = [0.5, 0.5, 0.5], scale = 0.2, size = (N, 3))
        vel = np.random.random(size = (N, 3)) * 10.0 - 5.0
        with h5py.File(fn, "w") as f:
            f.create_dataset("/particle_positions", data = pos)
            f.create_dataset("/particle_velocities", data = vel)
            f.create_dataset("/particle_masses", data = np.ones(N))
    return {'actions': [_generate],
            'targets': [fn],
            'uptodate': [run_once]}
