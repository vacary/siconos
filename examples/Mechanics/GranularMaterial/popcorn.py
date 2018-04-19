#!/usr/bin/env python

from siconos.io.mechanics_run import MechanicsHdf5Runner
import siconos.numerics as Numerics

from shutil import copyfile

fn = 'popcorn_scene_8_cubes.hdf5'
#fn = 'popcorn_scene_4_cubes.hdf5'
#fn = 'popcorn_scene_2_cubes.hdf5'

fn_simulation = fn.replace('scene', 'simulation')

print(fn, fn_simulation)

copyfile(fn, fn_simulation)


with MechanicsHdf5Runner(mode='r+', io_filename=fn_simulation, collision_margin =0.01) as io:
    io.run(t0=0,
           T=1,
           h=1e-4,
           multipoints_iterations=True,
           theta=1.0,
           Newton_max_iter=1,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=1000,
           tolerance=1e-3,
           exit_tolerance=10.0,
           output_frequency=10,
           violation_verbose=True,
           numerics_verbose=False)
