"""
AWP | Astrodynamics with Python by Alfonso Gonzalez
Illustration of pure rotation in velocity vector that changes angular momentum
    and eccentricity of an orbit.
    https://www.youtube.com/watch?v=rNpnzNKQrNg&list=PLOIRBaljOV8gn074rWFWYP1dCr2dJqWab

https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

3Blue1Brown's Summer of Math Exposition video entry
and Orbital Mechanics with Python 44
Gravity Assist Design via V-Infinity Matching to Explore
Our Solar System


Note: for select code, to prevent auto formatting (using vscode black),
    use the "# fmt: off" and "# fmt: on" commands.
"""
import sys  # needed to fix python importing issue

sys.path.append("src/python_tools")  # needed to fix python importing issue

# Libraries
import numerical_tools as nt
import numpy as np
import orbit_calculations as oc
import planetary_data as pd
import plotting_tools as pt
from Spacecraft import Spacecraft as SC

# fmt: off
if __name__ == "__main__":
    periapsis = pd.earth["radius"] + 4000.0  # km
    coes      = [periapsis / 0.3, 0.7, 0, 30.0, 0, 0]
    state     = oc.coes2state(coes)
    v_rot     = np.dot(nt.Cz(40.0 * nt.d2r), state[3:])
    state_rot = np.concatenate((state[:3], v_rot))

    sc0      = SC({"orbit_state": state, "tspan": "1"})
    sc_rot   = SC({"orbit_state": state_rot, "tspan": "1"})

    pt.plot_orbits(
        [sc0.states[:, :3], sc_rot.states[:, :3]],
        {
            "labels"      : ["Before", "After"],
            "colors"      : ["c", "m"],
            "azimuth"     : -90.0,
            "elevation"   : 90.0,
            "axes_custom" : 36000.0,
            "traj_lws"    : 2.5,
            "show"        : True,
        },
    )
