"""
AWP | Astrodynamics with Python by Alfonso Gonzalez
Earth-Venus-Mars-Earth gravity assist interplanetary trajectory implementation.
Orbital Mechanics with Python 44

https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

3Blue1Brown's Summer of Math Exposition video entry and,
Gravity Assist Design via V-Infinity Matching to Explore Our Solar System.

Earth-Venus-Mars-Earth gravity solved by the Interplanetary Trajectory
V-Infinity Matcher (ITVIM); called by other scripts in this directory.

Note: for select code, to prevent auto formatting (using vscode black),
    use the "# fmt: off" and "# fmt: on" commands.
"""

import sys  # needed to fix python importing issue

sys.path.append("src/python_tools")  # needed to fix python importing issue
# Libraries
import planetary_data as pd
import spice_data as sd
import spiceypy as spice
from ITVIM import ITVIM


def calc_EVME_1963():
    """
    Earth-Venus-Mars-Earth (EVME) 2 year trajectory, launched in 1966.
    Example from Richard Battin "Astronautical Guidance".
    """
    spice.furnsh(sd.leapseconds_kernel)
    spice.furnsh(sd.de432)

    sequence = [
        {"planet": 3, "time": "1966-02-10", "tm": -1},
        {
            "planet": 2,
            "planet_mu": pd.venus["mu"],
            "time": "1966-07-07",
            "tm": 1,
            "tol": 1e-5,
        },
        {
            "planet": 4,
            "planet_mu": pd.mars["mu"],
            "time": "1967-01-10",
            "tm": -1,
            "tol": 1e-5,
        },
        {"planet": 3, "planet_mu": pd.earth["mu"], "time": "1967-12-18", "tol": 1e-5},
    ]
    itvim = ITVIM({"sequence": sequence})
    itvim.print_summary()

    return itvim


if __name__ == "__main__":
    calc_EVME_1963()
