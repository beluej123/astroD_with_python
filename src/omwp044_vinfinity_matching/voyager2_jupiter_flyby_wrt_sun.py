"""
AWP | Astrodynamics with Python by Alfonso Gonzalez
Gravity Assist Design via V-Infinity Matching.
Create velocity vectors vs. time plot of Voyager 2 flyby w.r.t. Sun.

https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

3Blue1Brown's Summer of Math Exposition video entry
and Orbital Mechanics with Python 44

*********** NOTE ************
Run this script from the top project directory to get correct paths.
    Or you can change the paths to fit your needs.

Note: for select code, to prevent auto formatting (using vscode black),
    use the "# fmt: off" and "# fmt: on" commands.
"""
import sys  # needed to fix python importing issue

sys.path.append("src/python_tools")  # needed to fix python importing issue
# Libraries
import numpy as np
import plotting_tools as pt
import spice_data as sd
import spice_tools as st
import spiceypy as spice
from numerical_tools import norm

if __name__ == "__main__":
    # fmt: off
    spice.furnsh(sd.leapseconds_kernel)
    spice.furnsh("src/omwp044_vinfinity_matching/voyager2_jupiter_flyby.bsp")

    et        = spice.str2et("1979-07-09 TDB")
    dt        = 20 * 24 * 3600.0
    ets       = np.arange(et - dt, et + dt, 5000.0)
    states    = st.calc_ephemeris(-32, ets, "ECLIPJ2000", 10)
    state_jup = spice.spkgeo(5, et, "ECLIPJ2000", 10)[0]
    hline     = {"val": norm(state_jup[3:]), "color": "C3"}

    pt.plot_velocities(
        ets,
        states[:, 3:],
        {
            "time_unit": "days",
            "hlines": [hline],
            "show": True,
        },
    )
