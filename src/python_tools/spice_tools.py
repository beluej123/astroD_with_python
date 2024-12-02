"""
AWP | Astrodynamics with Python by Alfonso Gonzalez
SPICE convenience functions using SpiceyPy.
    https://www.youtube.com/watch?v=oJ8bBrtwZLk
    https://naif.jpl.nasa.gov/naif/spiceconcept.html
    https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/SPICE_Tutorials_all.pdf
Horizons web-ephemeris, https://ssd.jpl.nasa.gov/horizons/app.html#/
    

https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

S- Spacecraft ephemeris, given as a function of time (SPK; *.bsp files).
P- Planet, satellite, comet, or asteroid ephemerides (SPK),
    P component also logically includes certain physical, dynamical
    and cartographic constants for target bodies (PCK).
I- Instrument information containing descriptive data peculiar to
    the geometric aspects of a particular scientific instrument (IK).
C- Orientation information, containing a transformation, traditionally
    called the "C-matrix," which provides time-tagged pointing (orientation)
    angles.
E- Events information, summarizing mission activities (EK)- both planned
    and unanticipated - rarely used.

dsk, Digital Shape Kernel (for modeling the shape of a few natural bodies).
fk, Frames Kernel (for specific missions, it's not clear what they contain).
lsk, Leapseconds Kernel (for accurate time calculations).
pck, Planetary Constants Kernel (data related to orbits of major natural solar system bodies).
spk, Spacecraft and Planet Kernel (orbits and other details of planets, natural satellites)

Note: for select code, to prevent auto formatting (using vscode black),
    use the "# fmt: off" and "# fmt: on" commands.
"""

# Libraries
import os

import numpy as np
import spiceypy as spice
from numpy import array, zeros


# fmt: on
def calc_ephemeris(target, ets, frame, observer):
    """
    Convenience wrapper for spkezr and spkgeo
    """

    if type(target) == str:
        return array(spice.spkezr(target, ets, frame, "NONE", observer)[0])

    else:
        n_states = len(ets)
        states = zeros((n_states, 6))
        for n in range(n_states):
            states[n] = spice.spkgeo(target, ets[n], frame, observer)[0]
        return states


# fmt: off
def write_bsp(ets, states, args={}):
    """
    Write or append to a BSP / SPK kernel from a NumPy array
    """
    _args = {
        "bsp_fn"   : "traj.bsp",
        "spice_id" : -999,
        "center"   : 399,
        "frame"    : "J2000",
        "degree"   : 5,
        "verbose"  : True,
        "new"      : True,
        "comments" : "",
    }
    for key in args.keys():
        _args[key] = args[key]

    if _args["new"]:
        handle = spice.spkopn(_args["bsp_fn"], "SPK_file", len(_args["comments"]))
        action = "Wrote"
    else:
        handle = spice.spkopa(_args["bsp_fn"])
        action = "Updated"

    spice.spkw09(
        handle,
        _args["spice_id"],
        _args["center"],
        _args["frame"],
        ets[0],
        ets[-1],
        "0",
        _args["degree"],
        len(ets),
        states.tolist(),
        ets.tolist(),
    )

    spice.spkcls(handle)

    if _args["verbose"]:
        print(f'{action} { _args[ "bsp_fn" ] }.')
    return None


# fmt: on
def get_objects(filename, display=False):
    """
    Get spice objects

    Input Parameters:
    ----------
        filename : includes path
        display  : bool

    Returns:
    ----------
        ids     :
        names   :
        tcs_sec : time coverage's
        tcs_cal : time coverage's with calendar days
    """
    objects = spice.spkobj(filename)
    ids, names, tcs_sec, tcs_cal = [], [], [], []
    n = 0
    if display:
        print(f"\nObjects in {filename}")
    for o in objects:
        ids.append(o)
        # time coverage in seconds since J2000
        #   wnfetd= fetch a particular interval from a double precision window
        tc_sec = spice.wnfetd(spice.spkcov(filename, ids[n]), n)
        # convert time coverage to human readable
        tc_cal = [
            spice.timout(f, "YYYY-MON-DD HR:MN:SC.### (TDB) ::TDB") for f in tc_sec
        ]

        # append time coverages to putput lists
        tcs_sec.append(tc_sec)  # time coverage in seconds
        tcs_cal.append(tc_cal)
        # get body name
        try:
            names.append(id2body(o))
        except:
            # if body name does not exist
            names.append("Unknown Name")
        if display:
            # remember, time coverage start then time coverage finish
            print(
                f"id:{ids[-1]} \t{names[-1]}, \tvalid times: {tc_cal[0]} --> {tc_cal[1]}"
            )
    return ids, names, tcs_sec, tcs_cal

# create time array for time coverages
def tc2array(tcs, steps):
    arr=np.zeros((steps,1))
    arr[:,0] = np.linspace(tcs[0], tcs[1], steps)
    return arr


def id2body(id_):
    """function name is more intuitive than spice function name"""
    return spice.bodc2n(id_)


def find_directory(dir_name, start_dir=os.getcwd()):
    """
    Searches for a directory in parent directories.
    Input Parameters:
    ----------
        dir_name  : directory name you're looking for; must be parent, not subdirectory of parent
        start_dir : The directory to start searching from (defaults to current directory).

    Returns:
    ----------
        The absolute path to the directory if found, otherwise None.
    """
    current_dir = start_dir

    while True:
        if dir_name in os.listdir(current_dir):
            return os.path.join(current_dir, dir_name)

        parent_dir = os.path.dirname(current_dir)
        if parent_dir == current_dir:  # Reached the root directory
            return None

        current_dir = parent_dir

    return None


def test_spice_1():
    # solar system orbits, https://www.youtube.com/watch?v=NhnowBBLtmo
    # show spiceypy version
    print(f"{spice.tkvrsn('TOOLKIT')}")

    # find path to spice_data; search thru parents
    #   spice_data is a directory under OrbMech
    base_dir = find_directory(dir_name="OrbMech")
    if base_dir:
        print(f"Directory found: {base_dir}")
        # under the base directory is consistant structure holding spice data
        dir1 = os.path.join(base_dir, "!Orb_data", "spice_data", "naif0012.tls")
        dir2 = os.path.join(base_dir, "!Orb_data", "spice_data", "de421.bsp")
        print(f"dir1= {dir1}")
    else:
        print(f"Directory not found.")

    spice.furnsh(dir1)
    spice.furnsh(dir2)

    # Define the time of interest
    et = spice.str2et("2024-11-29T00:00:00")
    # Get the position of Earth relative to Sun
    # position, light_time = spice.spkpos('EARTH', et, 'J2000', 'NONE', 'SUN')
    position, light_time = spice.spkpos(
        targ="EARTH", et=et, ref="ECLIPJ2000", abcorr="NONE", obs="SUN"
    )

    print("Earth position relative to Sun:", position)

    # beginning coverage for Voyager 2 SPK kernel
    # et0 = spice.str2et("1977 AUG 20 15:32:32.182")


# fmt: on
if __name__ == "__main__":
    """
    To load SPICE kernels from any directory, set the os environment path.
        The original code author set an environment variable named 'AWP' to be
        the absolute path to the SPICE kernels from any directory.  Otherwise
        run this script from the src/directory.
        For example: AWP=/home/alfonso/pub/AWP
    """
    test_spice_1()  # explore spice files
