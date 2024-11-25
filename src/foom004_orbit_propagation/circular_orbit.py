"""
AWP | Astrodynamics with Python by Alfonso Gonzalez
Circular orbit propagation.

Fundamentals of Orbital Mechanics 4.
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Note: for select code, to prevent auto formatting (using vscode black),
    use the "# fmt: off" and "# fmt: on" commands.
"""

import sys  # needed to fix some python importing issues

sys.path.append("src/python_tools")  # fix some python importing issues

# 3rd party libraries
import numpy as np
import ode_tools as ot

# AWP library
import orbit_calculations as oc
import planetary_data as pd
import plotting_tools as pt

if __name__ == "__main__":
    # to help assignment alignments, temporarly turn off auto formatting
    # fmt: off
    r0_norm   = pd.earth["radius"] + 450.0  # km
    v0_norm   = (pd.earth["mu"] / r0_norm) ** 0.5  # km / s
    statei    = [r0_norm, 0, 0, 0, v0_norm, 0]
    tspan     = 100.0 * 60.0  # seconds
    dt        = 100.0  # seconds
    steps     = int(tspan / dt)
    ets       = np.zeros((steps, 1))
    states    = np.zeros((steps, 6))
    states[0] = statei

    # fmt: on
    for step in range(steps - 1):
        states[step + 1] = ot.rk4_step(oc.two_body_ode, ets[step], states[step], dt)

    pt.plot_orbits([states], {"show": True})
