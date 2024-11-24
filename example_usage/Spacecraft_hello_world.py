"""
AWP | Astrodynamics with Python by Alfonso Gonzalez
Two-body propagation with J2 perturbation for 100 periods.

https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering
"""

import sys

# sys.path() helps python find files, but does not fix the linting problem!
sys.path.append("src/python_tools")  # needed to fix python importing issue


from planetary_data import earth
from Spacecraft import Spacecraft as SC

if __name__ == "__main__":
    # coes= classic orbital elements
    #   orbit radius, eccentricity, inclination, true anomaly,
    coes = [earth["radius"] + 1000, 0.05, 30.0, 0.0, 0.0, 0.0]

    # pass in a dictionary defining calculations type
    sc = SC(
        {
            "coes"       : coes,  # sc initial condition
            "tspan"      : "100",  # str= no. of periods; or #=seconds
            "dt"         : 100.0,  # time steps in seconds
            "orbit_perts": {"J2": True},  # orbit pertabations added
            "propagate"  : True
        }
    )
    # a use-case to wait to propagate; in sc instance set 'propagate: False'
    #   then run seperate sc.propagate(), shown below
    # sc.propagate_orbit()
    
    # plots will show up in seperate "windows"
    sc.plot_3d({"show": True})
    sc.plot_coes()
    # sc.plot_states() # works fine, not useful for this mission
    # sc.plot_groundtracks() # 2024-11-23, not working
