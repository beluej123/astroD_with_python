"""
AWP | Astrodynamics with Python by Alfonso Gonzalez
Gravity Assist Design via V-Infinity Matching.
Voyager 2, 3D plot and orbital states plot with respect to Jupiter.

https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

3Blue1Brown's Summer of Math Exposition video;
    Orbital Mechanics with Python 44.

*********** NOTE ************
Run this script from the top project directory to get correct paths.
    Or you can change the paths to fit your needs.

Note: for select code, to prevent auto formatting (using vscode black),
    use the "# fmt: off" and "# fmt: on" commands.
"""
import sys  # needed to fix python importing issue

sys.path.append("src/python_tools")  # needed to fix python importing issue

#Library's
import numpy as np
import planetary_data as pd
import plotting_tools as pt
import spice_data as sd
import spice_tools as st
import spiceypy as spice
from numerical_tools import norm

# fmt:off
if __name__ == "__main__":
    spice.furnsh(sd.leapseconds_kernel)
    spice.furnsh("src/omwp044_vinfinity_matching/voyager2_jupiter_flyby.bsp")

    et     = spice.str2et("1979-07-09 TDB")
    dt     = 20 * 24 * 3600.0
    ets    = np.arange(et - dt, et + dt, 5000.0)
    states = st.calc_ephemeris(-32, ets, "ECLIPJ2000", 5)

    pt.plot_orbits(
        [states[:, :3]],
        {
            "labels"    : ["Voyager 2"],
            "colors"    : ["m"],
            "cb_radius" : pd.jupiter["radius"],
            "cb_cmap"   : "Oranges",
            "dist_unit" : "JR",
            "azimuth"   : -57,
            "elevation" : 15,
            "axes_mag"  : 0.5,
            "show"      : True,
        },
    )

    hline = {"val": norm(states[0, 3:]), "color": "m"}
    pt.plot_velocities(
        ets, states[:, 3:], {"time_unit": "hours", "hlines": [hline], "show": True}
    )
