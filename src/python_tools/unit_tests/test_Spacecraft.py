"""
AWP | Astrodynamics with Python by Alfonso Gonzalez
Spacecraft Class Unit Tests.

https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Note: for select code, to prevent auto formatting (using vscode black),
    use the "# fmt: off" and "# fmt: on" commands.
"""

import sys  # needed to fix python importing issue

sys.path.append("src/python_tools")  # needed to fix python importing issue
# Libraries
import numpy as np
import planetary_data as pd
import plotting_tools as pt
import pytest
import spice_data as sd
import spice_tools as st
import spiceypy as spice
from Spacecraft import Spacecraft as SC

# Treat all warnings as errors
pytestmark = pytest.mark.filterwarnings("error")

# fmt: off
def test_Spacecraft_basic_propagation(plot=False):
    sc = SC(
        {
            "coes": [pd.earth["radius"] + 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            "tspan": "1",
            "dt": 100.0,
            "rtol": 1e-8,
        }
    )
    assert sc.cb == pd.earth

    """
	Since inclination is 0, all z-axis components of Spacecraft
	position and velocity should be 0
	"""
    assert np.all(sc.states[:, 2] == 0.0)
    assert np.all(sc.states[:, 5] == 0.0)

    """
	Since there are no orbital perturbations, all COEs except
	true anomaly should	be close to constant
	(not exactly constant due to numerical error).
	However, since this is a circular orbit, periapsis is loosely
	defined, causing errors in true anomaly (not exactly linear),
	argument of periapsis (bounces back and forth between
	0 and 359.9 degrees), eccentricity (random noise), and
	semi-major axis (drift / random noise)
	"""
    sc.calc_coes()
    assert sc.coes_calculated
    assert pytest.approx(sc.coes[:, 0], abs=1e-3) == pd.earth["radius"] + 1000.0  # sma
    assert pytest.approx(sc.coes[:, 1], abs=1e-6) == 0.0  # ecc
    assert np.all(sc.coes[:, 2]                   == 0.0)  # inc
    assert np.all(sc.coes[:, 5]                   == 0.0)  # raan

    sc.calc_apoapses_periapses()
    apse_diffs = sc.apoapses - sc.periapses
    assert pytest.approx(apse_diffs, abs=1e-3) == 0.0

    if plot:
        sc.plot_coes()


def test_Spacecraft_inclination_latitude(plot=False):
    spice.furnsh(sd.pck00010)

    sc = SC(
        {
            "coes": [pd.earth["radius"] + 1000.0, 0.01, 50.0, 0.0, 0.0, 0.0],
            "tspan": "3",
            "dt": 10.0,
            "rtol": 1e-8,
        }
    )
    assert sc.cb == pd.earth

    """
	Given that this spacecraft has 50 degrees inclination,
	the latitude coordinates should remain in between
	-50 and 50 degrees, since inclination is defined as
	the angle between the orbital plane and Earth's
	equatorial plane
	"""
    sc.calc_latlons()
    assert np.all(sc.latlons[:, 2] <= 50.0) and np.all(sc.latlons[:, 2] >= -50.0)

    if plot:
        sc.plot_groundtracks()


def test_Spacecraft_minimum_altitude_stop_condition(plot=False):
    sc = SC(
        {
            "coes"            : [pd.earth["radius"] + 1000.0, 0.5, 0.0, 90.0, 0.0, 0.0],
            "tspan"           : "1",
            "dt"              : 100.0,
            "rtol"            : 1e-9,
            "stop_conditions" : {"min_alt": 100.0},  # km
        }
    )
    assert sc.cb == pd.earth

    """
	Since inclination is 0, all z-axis components of Spacecraft
	position and velocity should be 0
	"""
    assert np.all(sc.states[:, 2] == 0.0)
    assert np.all(sc.states[:, 5] == 0.0)

    """
	The orbital elements were set such that minimum altitude
	crossing would occur, triggering the min_alt stop condition
	"""
    assert sc.ode_sol.success
    assert sc.ode_sol.message == "A termination event occurred."
    assert sc.ode_sol.status  == 1

    assert (
        pytest.approx(
            np.linalg.norm(sc.ode_sol.y_events[0][:3]) - sc.cb["radius"] - 100.0,
            abs=6.0e-3,
        )
        == 0.0
    )

    if plot:
        sc.plot_altitudes(
            {
                "hlines"    : [{"val": 100.0, "color": "c"}],
                "time_unit" : "hours",
                "title"     : "test_Spacecraft_minimum_altitude_stop_condition",
                "show"      : True,
            }
        )


def test_Spacecraft_maximum_altitude_stop_condition(plot=False):
    stop_conditions = {"max_alt": pd.earth["SOI"] - pd.earth["radius"]}
    sc = SC(
        {
            "coes"            : [-(pd.earth["radius"] + 1000.0), 1.8, 0.0, 10.0, 0.0, 0.0],
            "tspan"           : 5 * 24 * 3600.0,
            "rtol"            : 1e-9,
            "stop_conditions" : stop_conditions,
        }
    )
    assert sc.cb == pd.earth

    """
	Since inclination is 0, all z-axis components of Spacecraft
	position and velocity should be 0
	"""
    assert np.all(sc.states[:, 2] == 0.0)
    assert np.all(sc.states[:, 5] == 0.0)

    """
	The orbital elements were set such that maximum altitude
	crossing would occur, triggering the max_alt stop condition
	"""
    assert sc.ode_sol.success
    assert sc.ode_sol.message == "A termination event occurred."
    assert sc.ode_sol.status  == 1

    assert (
        pytest.approx(
            np.linalg.norm(sc.ode_sol.y_events[1][:3]) - sc.cb["radius"], abs=1e-3
        )
        == sc.cb["SOI"] - sc.cb["radius"]
    )

    if plot:
        hline = {"val": pd.earth["SOI"] - pd.earth["radius"], "color": "c"}
        sc.plot_altitudes(
            {
                "hlines"    : [hline],
                "time_unit" : "hours",
                "title"     : "test_Spacecraft_maximum_altitude_stop_condition",
                "show"      : True,
            }
        )


def test_Spacecraft_enter_SOI_voyager2_jupiter(plot=False):
    """
    Horizons web-ephemeris, https://ssd.jpl.nasa.gov/horizons/app.html#/
        Earth: Ecliptic: Solar System Barycenter
        2443376.147594699 = A.D. 1977-Aug-20 15:32:32.1820 TDB 
        X = 1.284144508572808E+08 Y =-8.128705432252860E+07 Z =-1.253680534372479E+04
        VX= 1.539202448857433E+01 VY= 2.510752773341791E+01 VZ=-1.147795727494128E-04
        LT= 5.069500111576398E+02 RG= 1.519797899280763E+08 RR=-4.234681274407071E-01
    """
    spice.furnsh(sd.leapseconds_kernel)
    spice.furnsh(sd.de432)
    frame = "ECLIPJ2000"

    # beginning of coverage for Voyager 2 SPK kernel
    et0 = spice.str2et("1977 AUG 20 15:32:32.182")

    """
	Voyager 2 initial state w.r.t Earth at et0
	"""
    state0 = [
        7.44304239e03, -5.12480267e02, 2.37748995e03,  # [km]
        5.99287925e00,  1.19056063e01, 5.17252832e00,  # [km/s]
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
            "tspan"          : 100000, # [s]
            "stop_conditions": stop_conditions,
        }
    )
    assert sc0.cb              == pd.earth
    assert sc0.ode_sol.success
    assert sc0.ode_sol.message == "A termination event occurred."
    assert sc0.ode_sol.status  == 1

    """
	Calculate spacecraft state w.r.t solar system barycenter
	at ephemeris time when spacecraft left Earth SOI
	"""
    state_earth = spice.spkgeo(399, sc0.ets[-1], frame, 0)[0]
    state1      = sc0.states[-1, :6] + state_earth

    """
	Now model the spacecraft as a heliocentric elliptical orbit
	and propagate until enter Jupiter SOI
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
    assert sc1.cb              == pd.sun
    assert sc1.ode_sol.success
    assert sc1.ode_sol.message == "A termination event occurred."
    assert sc1.ode_sol.status  == 1

    if plot:
        ets        = np.concatenate((sc0.ets, sc1.ets))
        rs_earth   = st.calc_ephemeris(3, ets, frame, 10)[:, :3]
        rs_jupiter = st.calc_ephemeris(5, ets, frame, 10)[:, :3]
        rs0        = sc0.states[:, :3] + st.calc_ephemeris(3, sc0.ets, frame, 10)[:, :3]

        labels = ["$Voyager2_{EarthSOI}$", "$Voyager2_{SunSOI}$", "Earth", "Jupiter"]

        pt.plot_orbits(
            [rs0, sc1.states[:, :3], rs_earth, rs_jupiter],
            {"labels": labels, "colors": ["m", "c", "b", "C3"], "show": True},
        )


if __name__ == "__main__":
    test_Spacecraft_basic_propagation(plot=True)
    # test_Spacecraft_inclination_latitude( plot = True ) # fix cities problem
    test_Spacecraft_minimum_altitude_stop_condition(plot=True)
    test_Spacecraft_maximum_altitude_stop_condition(plot=True)
    test_Spacecraft_enter_SOI_voyager2_jupiter(plot=True)
