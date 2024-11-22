"""
AWP | Astrodynamics with Python inspired by Alfonso Gonzalez
Hohmann transfer calculations.
https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering
Hohmann: Python 36 & Python 28
https://www.youtube.com/watch?v=p4fVnCbZxf0&list=PLOIRBaljOV8gn074rWFWYP1dCr2dJqWab&index=14
https://www.youtube.com/watch?v=35eQ9FHom7o

Orbital Mechanics began with video Python 28
2024-11-20+ Jeff Belue creation.
For code readability, prevent some specific auto formatting, I use the
    "# fmt: off" and "# fmt: on" commands.  I use black formatter within vscode.

References:
----------
    See references.py for references list.
"""

import sys  # needed to fix python importing issue

sys.path.append("src/python_tools")  # needed to fix python importing issue
import orbit_calculations as orbit_c
import planetary_data as pd
import plotting_tools as pt
import Spacecraft as SC


def test_hohmann_1():
    cb = pd.earth  # central body
    
    # initial and final altitudes
    r0 = 1000  # [km]
    r1 = 10000  # [km]

    # coe (classical orbital elements) of initial and final orbits
    #   coe: ecc, incl, TA, argP, raan
    # equatorial case
    coes0 = [cb["radius"] + r0, 0.0, 0.0, 0.0, 0.0, 0.0]
    coes1 = [cb["radius"] + r1, 0.0, 0.0, 180.0, 0.0, 0.0]
    # inclined/rotated case
    coes0_rot = [cb["radius"] + r0, 0.0, 30.0, 0.0, 100.0, 20.0]
    coes1_rot = [cb["radius"] + r1, 0.0, 30.0, 180.0, 100.0, 20.0]
    # call scalar hohmann transfer function
    delta_vs_scalar, t_transfer_scalar = orbit_c.hohmann_transfer_scalars(r0, r1)

    # print hohmann transfer scalar parameters to terminal
    print(f"delta V0:\t{delta_vs_scalar[0]:.5g} [km/s]")
    print(f"delta V1:\t{delta_vs_scalar[1]:.5g} [km/s]")
    print(
        f"transfer time:\t{t_transfer_scalar:.5g} [s], {t_transfer_scalar/3600:.5g} [hr]"
    )

    # **** troubleshooting, below **********
    sc_transfer, delta_vs = orbit_c.hohmann_transfer(
        coes0=coes0, coes1=coes1, propagate=True
    )
    print(f"sc_transfer= {sc_transfer}")
    # **** troubleshooting, above **********
    
    # sc0, sc1, sc_transfer, delta_vs = orbit_c.hohmann_transfer(
    #     coes0=coes0, coes1=coes1, propagate=True
    # )
    # sc0_rot, sc1_rot, sc_transfer_rot, delta_vs_rot = orbit_c.hohmann_transfer(
    #     coes0=coes0, coes1=coes1, propagate=True
    # )

    # pt.plot_orbits([sc0.rs, sc1.rs, sc_transfer.rs, sc0_rot.rs, sc1_rot.rs, sc_transfer_rot.rs],
    #               labels=['initial', 'final', 'transfer', 'initial_', 'final_', 'transfer_'],
    #               az=-45.0, el=0.0, axes=9500,
    #               save_plot=True, no_axes=True, title='', dpi=300,
    #               output_dir='v27_28_hohmann'):


if __name__ == "__main__":
    test_hohmann_1()
