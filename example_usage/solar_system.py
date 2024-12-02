"""
AWP | Astrodynamics with Python by Alfonso Gonzalez
Plot solar system planets orbits; check video time at 21:38+
    https://www.youtube.com/watch?v=NhnowBBLtmo

https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering
"""

import sys

# sys.path() helps python find files, but does not fix the linting problem!
sys.path.append("src/python_tools")  # needed to fix python importing issue
import os

import matplotlib.pyplot as plt
import numpy as np
import spiceypy as spice

plt.style.use("dark_background")

import planetary_data as pd
import plotting_tools as pt
import spice_tools as st
from Spacecraft import Spacecraft as SC  # orbit propagator is in here

# central body
cb = pd.earth
# total steps for ephemeris data
STEPS = 100000
# planets reference frame & planetary observer
FRAME = "ECLIPJ2000"
OBSERVER = "SUN"

if __name__ == "__main__":
    # find path to spice_data; search thru parents
    #   spice_data is a directory under OrbMech
    base_dir = st.find_directory(dir_name="OrbMech")
    if base_dir:
        print(f"Parent directory found: {base_dir}")
        # under the base directory is consistant structure holding spice data
        dir1 = os.path.join(base_dir, "!Orb_data", "spice_data", "naif0012.tls")
        dir2 = os.path.join(base_dir, "!Orb_data", "spice_data", "de432s.bsp")
        # print(f"dir1= {dir1}")
    else:
        print(f"Directory not found.")

    # TODO make spice meta data for solar system ephemeris
    spice.furnsh(dir1)
    spice.furnsh(dir2)

    ids, names, tcs_sec, tcs_cal = st.get_objects(filename=dir2, display=True)
    # only include barycenters
    names = [f for f in names if "BARYCENTER" in f]

    # create time array for ephemeris data for all bodies
    times = st.tc2array(tcs_sec[0], STEPS)

    # empty data list
    rs = []
    # for each body in the solar system
    for name in names:
        # add ephemeris data to list
        # rs.append(st.get_ephemeris_data(name, times, FRAME, OBSERVER))
        rs.append(st.calc_ephemeris(name, times, FRAME, OBSERVER))

    # check arguements: plot_orbits(rs, args, vectors=[]):
    pt.plot_orbits(
        rs,
        {
            "traj_lws": 1,
            "show": True,
        },
    )
