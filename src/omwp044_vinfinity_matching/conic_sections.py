"""
AWP | Astrodynamics with Python by Alfonso Gonzalez
Plot the conic sections; circle, ellipse, parabola, hyperbola.

https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering
https://www.youtube.com/@alfonsogonzalez-astrodynam2207/featured

3Blue1Brown's Summer of Math Exposition video entry
and Orbital Mechanics with Python 44
Gravity Assist Design via V-Infinity Matching to Explore Our Solar System

Note: for select code, to prevent auto formatting (using vscode black),
    use the "# fmt: off" and "# fmt: on" commands.
"""

import sys  # needed to fix python importing issue

sys.path.append("src/python_tools")  # needed to fix python importing issue
# Libraries
import orbit_calculations as oc
import planetary_data as pd
import plotting_tools as pt
import spiceypy as spice
from Spacecraft import Spacecraft as SC

# fmt: off
if __name__ == "__main__":
    periapsis = pd.earth["radius"] + 4000.0  # km

    coes_circular   = [periapsis,    0,   0, 0, 0, 0]
    coes_elliptical = [periapsis / 0.3, 0.7, 0, 0, 0, 0]
    state_parabolic = spice.conics(
        [periapsis, 1.0, 0, 0, 0, -10.0, 0, pd.earth["mu"]], 0
    )
    state_hyperbolic = spice.conics(
        [periapsis, 2.5, 0, 0, 0, -10.0, 0, pd.earth["mu"]], 0
    )

    sc_circular   = SC({"coes": coes_circular, "tspan": "1"})
    sc_elliptical = SC({"coes": coes_elliptical, "tspan": "2"})
    sc_parabolic  = SC({"orbit_state": state_parabolic, "tspan": 70000.0})
    sc_hyperbolic = SC({"orbit_state": state_hyperbolic, "tspan": 30000.0})

    rs = [
        sc_circular.states[:, :3],
        sc_elliptical.states[:, :3],
        sc_parabolic.states[:, :3],
        sc_hyperbolic.states[:, :3],
    ]
    labels = [
        "Circular $(ecc=0.0)$",
        "Elliptical $(ecc=0.7)$",
        "Parabolic $(ecc=1.0)$",
        "Hyperbolic $(ecc=3.0)$",
    ]

    sc_elliptical.plot_states({"time_unit": "hours", "show": True})

    pt.plot_orbits(
        rs,
        {
            "labels"     : labels,
            "colors"     : ["r", "g", "b", "m"],
            "azimuth"    : -90.0,
            "elevation"  : 90.0,
            "axes_custom": 36000.0,
            "hide_axes"  : True,
            "traj_lws"   : 2.5,
            "show"       : True,
        },
    )
