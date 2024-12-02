'''
AWP | Astrodynamics with Python by Alfonso Gonzalez
Planetary Data Library.

https://github.com/alfonsogonzalez/AWP
https://www.youtube.com/c/AlfonsoGonzalezSpaceEngineering

Orbital Elements Naming Collection:
Start with Kepler coe (classic orbital elements).
    https://ssd.jpl.nasa.gov/planets/approx_pos.html
    Horizons web-ephemeris, https://ssd.jpl.nasa.gov/horizons/app.html#/

    o_type : int  , [-] orbit type (python dictionary list)
                    0:"circular", 1:"circular inclined", 2:"circular equatorial"
                    3:"elliptical", 4:"elliptical equatorial"
                    5:"parabolic", 6:"parabolic equatorial"
                    7:"hyperbolic", 8:"hyperbolic equatorial"
    sp     : float, [km or au] semi-parameter (aka p)
    sma    : float, [km or au] semi-major axis (aka a)
    ecc    : float, [--] eccentricity
    incl   : float, [rad] inclination
    raan   : float, [rad] right ascension of ascending node,
                    also called Longitude of Ascending Node (Omega, or capital W)
    w_     : float, [rad] argument of periapsis (aka aop, or arg_p)
    TA     : float, [rad] true angle/anomaly (aka t_anom, or theta)

    alternative coe's including circular & equatorial:
    Lt0    : float, [rad] true longitude at epoch, circular equatorial
                    Position on the ecliptic, accounting for its inclination.
                    when incl=0, ecc=0
    w_bar  : float, [rad] longitude of periapsis (aka II), equatorial
                NOTE ** NOT argument of periapsis, w_ ??????????????????? **
                Note, w_bar = w_ + raan, measured in 2 planes (Vallado [4] p.1015)
    u_     : float, [rad] argument of lattitude (aka ), circular inclined

    Other orbital elements:
    w_p    : float [rad] longitude of periapsis (aka w_bar) ??
    L_     : float, [deg] mean longitude
                NOT mean anomaly, M
                L_ = w_bar + M
    wt_bar : float, [rad] true longitude of periapsis
                measured in one plane
    M_     : mean anomaly, often replaces TA
    t_peri : float, [jd] time of periapsis passage

    circular, e=0: w_ and TA = undefined;
        use argument of latitude, u_; u_=acos((n_vec X r_vec)/(n_mag * r_mag))
        If r_vec[2] < 0 then 180 < u < 360 degree

    equatorial, i=0 or 180 [deg]: raan and w_ = undefined
        use longitude of periapsis, II (aka w_bar); II=acos(e_vec[0]/e_mag)
        If e_vec[1] < 0 then 180 < II < 360 degree

    circular & equatorial, e=0 and i=0 or i=180: w_ and raan and TA = undefined;
        use true longitude, Lt0 = angle between r0 & I-axis; Lt0=acos(r_vec[1]/r_mag)
        If r_mag[1] < 0 then 180 < Lt0 < 360 degree

From JPL Horizizons, osculating elements:
    Symbol & meaning [1 au= 149597870.700 km, 1 day= 86400.0 s]:
    JDTDB  Julian Day Number, Barycentric Dynamical Time
    EC     Eccentricity, e
    QR     Periapsis distance, q (au)
    IN     Inclination w.r.t X-Y plane, i (degrees)
    OM     Longitude of Ascending Node, OMEGA, (degrees)
    W      Argument of Perifocus, w (degrees)
    Tp     Time of periapsis (Julian Day Number)
    N      Mean motion, n (degrees/day)
    MA     Mean anomaly, M (degrees)
    TA     True anomaly, nu (degrees)
    A      Semi-major axis, a (au)
    AD     Apoapsis distance (au)
    PR     Sidereal orbit period (day)
'''

# gravitational constant
G_meters = 6.67430e-11       # m**3 / kg / s**2
G        = G_meters * 10**-9 # km**3/ kg / s**2

venus = {
		'name'            : 'Venus',
		'spice_name'      : 'VENUS BARYCENTER',
		'SPICE_ID'        : 2,
		'mass'            : 4.867e24,
		'mu'              : 3.2485859200000006E+05,
		'radius'          : 6051.8,
		'sma'             : 108.209e6,   # km
		'SOI'             : 617183.2511, # km
		'deorbit_altitude': 100.0,       # km
		'cmap'            : 'Wistia',
		'body_fixed_frame': 'IAU_VENUS',
		'traj_color'      : 'y'
		}

earth = {
		'name'            : 'Earth',
		'spice_name'      : 'EARTH',
		'SPICE_ID'        : 399,
		'mass'            : 5.972e24,
		'mu'              : 5.972e24 * G,
		'radius'          : 6378.0,
		'J2'              : 1.081874e-3,
		'sma'             : 149.596e6, # km
		'SOI'             : 926006.6608, # km
		'deorbit_altitude': 100.0, # km
		'cmap'            : 'Blues',
		'body_fixed_frame': 'IAU_EARTH',
		'traj_color'      : 'b'
		}

moon = {
		'name'            : 'Moon',
		'spice_name'      : 'MOON',
		'SPICE_ID'        : 301,
		'mass'            : 5.972e24,
		'mu'              : 5.972e24 * G,
		'radius'          : 1737.4,
		'J2'              : 1.081874e-3,
		'sma'             : 149.596e6, # km
		'SOI'             : 926006.6608, # km
		'deorbit_altitude': 100.0, # km
		'cmap'            : 'Blues',
		'body_fixed_frame': 'IAU_EARTH',
		'traj_color'      : 'b'
		}

mars = {
		'name'            : 'Mars',
		'spice_name'      : 'MARS BARYCENTER',
		'SPICE_ID'        : 4,
		'mass'            : 6.39e23,
		'mu'              : 4.282837362069909E+04,
		'radius'          : 3397.0,
		'sma'             : 227.923e6, # km
		'SOI'             : 0.578e6,   # km
		'deorbit_altitude': 50.0,      # km
		'cmap'            : 'Reds',
		'body_fixed_frame': 'IAU_MARS',
		'traj_color'      : 'r'
		}

jupiter = {
		'name'            : 'Jupiter',
		'spice_name'      : 'JUPITER BARYCENTER',
		'SPICE_ID'        : 5,
		'mass'            : 1.898e27,
		'mu'              : 1.26686e8,
		'radius'          : 71490.0,   # km
		'sma'             : 778.570e6, # km
		'deorbit_altitude': 1000.0,    # km
		'SOI'             : 48.2e6,    # km
		'body_fixed_frame': 'IAU_JUPITER',
		'traj_color'      : 'C3'
}

io = {
		'name'            : 'Io',
		'spice_name'      : 'Io',
		'SPICE_ID'        : 501,
		'mass'            : 1.898e27,
		'mu'              : 5.959916033410404E+03,
		'radius'          : 1821.6,   # km
		'deorbit_altitude': 10.0,    # km
		'traj_color'      : 'C1'
}

europa = {
		'name'            : 'Europa',
		'spice_name'      : 'Europa',
		'SPICE_ID'        : 502,
		'mu'              : 3.202738774922892E+03,
		'radius'          : 1560.8,   # km
		'deorbit_altitude': 10.0,    # km
		'traj_color'      : 'C2'
}

ganymede = {
		'name'            : 'Ganymede',
		'spice_name'      : 'Ganymede',
		'SPICE_ID'        : 503,
		'mu'              : 9.887834453334144E+03,
		'radius'          : 2631.2,   # km
		'deorbit_altitude': 100.0,    # km
		'traj_color'      : 'C3'
}

callisto = {
		'name'            : 'Callisto',
		'spice_name'      : 'Callisto',
		'SPICE_ID'        : 504,
		'mu'              : 7.179289361397270E+03,
		'radius'          : 2410.3,   # km
		'deorbit_altitude': 10.0,    # km
		'traj_color'      : 'C4'
}

saturn = {
	'name'            : 'Saturn',
	'spice_name'      : 'SATURN BARYCENTER',
	'SPICE_ID'        : 6,
	'mass'            : 568.34e24,
	'radius'          : 58232.0,
	'mu'              : 37.931e6,
	'sma'             : 1433.529e6,
	'deorbit_altitude': 1000.0,
	'SOI'             : 54890347.727,
	'traj_color'      : 'C2'
}

sun = {
	'name'            : 'Sun',
	'SPICE_ID'        : 10,
	'mass'            : 1.989e30,
	'mu'              : 1.3271244004193938E+11,
	'radius'          : 695510.0,
	'deorbit_altitude': 1.2 * 695510.0,
	'cmap'            :'gist_heat'
}

bodies = [
	venus, earth, moon, mars, 
	jupiter, io, europa, ganymede, callisto,
	saturn, sun ]

for body in bodies:
	body[ 'diameter' ] = body[ 'radius' ] * 2
