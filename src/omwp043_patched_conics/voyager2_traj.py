"""
AWP | Astrodynamics with Python by Alfonso Gonzalez
Patched Conics Propagation and Spacecraft Propagation Stop Conditions.
    Retrieve Voyager 2 initial state, patched conics propagation.
    Orbital Mechanics with Python 43
    
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

from Voyager_2.m05016u.merged.bsp which can be found at:
    https://naif.jpl.nasa.gov/pub/naif/VOYAGER/kernels/spk/

Note: run this script from the top project directory to get correct paths.
    Or you can change the paths to fit your needs.

Note: for select code, to prevent auto formatting (using vscode black),
    use the "# fmt: off" and "# fmt: on" commands.

Horizons web-ephemeris, https://ssd.jpl.nasa.gov/horizons/app.html#/
    Earth: Ecliptic: Solar System Barycenter
    2443376.147594699 = A.D. 1977-Aug-20 15:32:32.1820 TDB 
    X = 1.284144508572808E+08 Y =-8.128705432252860E+07 Z =-1.253680534372479E+04
    VX= 1.539202448857433E+01 VY= 2.510752773341791E+01 VZ=-1.147795727494128E-04
    LT= 5.069500111576398E+02 RG= 1.519797899280763E+08 RR=-4.234681274407071E-01
"""

# Python libraries
import sys  # needed to fix python importing issue

sys.path.append("src/python_tools")  # needed to fix python importing issue
import os

import planetary_data as pd
import plotting_tools as pt
import spice_tools as st
import spiceypy as spice
from numpy import concatenate
from Spacecraft import Spacecraft as SC

frame = "ECLIPJ2000"
# fmt: off
if __name__ == "__main__":
    """
    To load SPICE kernels from any directory, set the os environment path.
        The original code author set an environment variable named 'AWP' to be
        the absolute path to the SPICE kernels from any directory.  Otherwise
        run this script from the src/ directory.
        For example: AWP=/home/alfonso/pub/AWP
    """
    AWP_path = os.environ.get("AWP")
    if AWP_path is not None:
        base_dir = AWP_path
    else:
        base_dir = "." # remember to start from program root !

    spice.furnsh(os.path.join(base_dir, "data/spice/lsk/naif0012.tls"))
    spice.furnsh(os.path.join(base_dir, "data/spice/spk/de432s-1977-1983.bsp"))

    # beginning coverage for Voyager 2 SPK kernel
    et0 = spice.str2et("1977 AUG 20 15:32:32.182")

    """
	Voyager 2 initial state w.r.t Earth at et0
	"""
    state0 = [
        7.44304239e03, -5.12480267e02, 2.37748995e03,  # [km]
        5.99287925e00, 1.19056063e01,  5.17252832e00,  # [km/s]
    ]

    """
	Propagate Voyager 2 trajectory until it exits Earth's sphere of influence,
	where maximum altitude stop condition will be triggered
	"""
    stop_conditions = {"max_alt": pd.earth["SOI"] - pd.earth["radius"]}
    sc0 = SC(
        {
            "orbit_state"    : state0, # r_vec, v_vec [km, km/s]
            "et0"            : et0,    # ephemeris time
            "frame"          : frame,
            "tspan"          : 100000, #[s]
            "dt"             : 10000,  #[s]
            "stop_conditions": stop_conditions,
        }
    )
    """
	Calculate spacecraft state w.r.t solar system barycenter
	    at ephemeris time when spacecraft left Earth SOI.
	"""
    state_earth = spice.spkgeo(399, sc0.ets[-1], frame, 0)[0]
    state1      = sc0.states[-1, :6] + state_earth
    """
	Now model the spacecraft as a heliocentric elliptical orbit
	    and propagate until enter Jupiter SOI.
	"""
    sc1 = SC(
        {
            "orbit_state"    : state1,
            "et0"            : sc0.ets[-1],
            "frame"          : frame,
            "tspan"          : 5 * 365 * 24 * 3600.0,
            "dt"             : 30000,
            "stop_conditions": {"enter_SOI": pd.jupiter},
            "cb"             : pd.sun,
        }
    )

    """
	Calculate spacecraft state w.r.t Jupiter
	    at ephemeris time when spacecraft reaches Jupiter SOI.
	"""
    state_jupiter = spice.spkgeo(5, sc1.ets[-1], frame, 0)[0]
    state2        = sc1.states[-1, :6] - state_jupiter
    """
	Now model the spacecraft as a Jovicentric hyperbolic orbit
	    and propagate until exit Jupiter SOI.
	"""
    stop_conditions = {"max_alt": pd.jupiter["SOI"] - pd.jupiter["radius"]}
    sc2 = SC(
        {
            "orbit_state"    : state2,
            "et0"            : sc1.ets[-1],
            "frame"          : frame,
            "tspan"          : 200 * 24 * 3600.0,
            "dt"             : 20000,
            "stop_conditions": stop_conditions,
            "cb"             : pd.jupiter,
        }
    )

    """
	Calculate spacecraft state w.r.t solar system barycenter
	    at ephemeris time when spacecraft left Jupiter SOI.
	"""
    state_jupiter = spice.spkgeo(5, sc2.ets[-1], frame, 0)[0]
    state3        = sc2.states[-1, :6] + state_jupiter
    
    """
	Now model the spacecraft as a heliocentric elliptical orbit
	    and propagate until enter Saturn SOI.
	"""
    sc3 = SC(
        {
            "orbit_state"    : state3,
            "et0"            : sc2.ets[-1],
            "frame"          : frame,
            "tspan"          : 3 * 365 * 24 * 3600.0,
            "dt"             : 30000,
            "stop_conditions": {"enter_SOI": pd.saturn},
            "cb"             : pd.sun,
        }
    )

    """
	Create a NumPy array with all ets (ephemeris times)
	"""
    ets = concatenate((sc0.ets, sc1.ets, sc2.ets, sc3.ets))

    states_earth   = st.calc_ephemeris(399, ets, frame, 0)[:, :3]
    states_jupiter = st.calc_ephemeris(5,   ets, frame, 0)[:, :3]
    states_saturn  = st.calc_ephemeris(6,   ets, frame, 0)[:, :3]
    labels = [
        "Earth-Centered",
        "Heliocentric",
        "Jovicentric",
        "Heliocentric",
        "Earth",
        "Jupiter",
        "Saturn",
    ]
    colors = ["m", "c", "m", "c", "b", "C3", "C1"]

    # print(f"states_saturn= {states_saturn}") # debug cmd
    # print(f"sc0.__dict__)= {sc0.__dict__}") # debug cmd
    # exit() # debugging tool; i cannot get vscode degugger to work as I need!

    # ensure all states are heliocentric
    rs0 = sc0.states[:, :3] + states_earth[: sc0.n_steps]
    rs1 = sc1.states[:, :3]
    c0  = sc0.n_steps + sc1.n_steps
    c1  = c0 + sc2.n_steps
    rs2 = sc2.states[:, :3] + states_jupiter[c0:c1]
    rs3 = sc3.states[:, :3]

    pt.plot_orbits(
        [rs0, rs1, rs2, rs3, states_earth, states_jupiter, states_saturn],
        {
            "labels"   : labels,
            "colors"   : colors,
            "axes_mag" : 1.0,
            "traj_lws" : 2,
            "azimuth"  : -90,
            "elevation": 94,
            "show"     : True,
        },
    )

    """
	The following part of this script was used to create the GIF
	showed in this video. At this time I'm not ready to show the
	OrbitalAnimator class since I am going to implement class inheritance
	to organize all the animator classes together, but I kept this
	part in here in case anyone is curious
	"""
    if False:
        from sys import path

        path.append("/home/alfonso/AWP/python_tools")
        from OrbitalAnimator import DummySC, animate_orbits

        count = sc0.step
        rs_1 = np.zeros((sc1.step + count, 3))
        for n in range(count):
            rs_1[n] = rs1[0]
        rs_1[count:] = rs1

        count += sc1.step
        rs_2 = np.zeros((sc2.step + count, 3))
        for n in range(count):
            rs_2[n] = rs2[0]
        rs_2[count:] = rs2

        count += sc2.step
        rs_3 = np.zeros((sc3.step + count, 3))
        for n in range(count):
            rs_3[n] = rs3[0]
        rs_3[count:] = rs3

        sc_0 = DummySC(ets, rs0)
        sc_1 = DummySC(ets, rs_1)
        sc_2 = DummySC(ets, rs_2)
        sc_3 = DummySC(ets, rs_3)
        sc_earth = DummySC(ets, states_earth)
        sc_jupiter = DummySC(ets, states_jupiter)
        sc_saturn = DummySC(ets, states_saturn)

        animate_orbits(
            [sc_0, sc_1, sc_2, sc_3, sc_earth, sc_jupiter, sc_saturn],
            {
                "cb_radius"   : pd.sun["radius"],
                "sc_labels"   : labels,
                "traj_colors" : colors,
                "axes_mag"    : 7,
                "freq"        : 50,
                "animation_fn": "/mnt/c/Users/alfon/AWP/omwp43_patched_conics/voyager2.gif",
            },
        )
