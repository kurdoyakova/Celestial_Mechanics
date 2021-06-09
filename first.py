from astropy import units as u
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
import numpy as np
from astropy.time import Time
from astropy.coordinates import ICRS

t = Time('2019-05-27')
t.format = 'jyear_str'

a = 1.769562 * u.AU
ecc = 0.087686 * u.one
inc = 3.31053 * u.deg
raan = 8.92533 * u.deg
argp = 320.11310 * u.deg
nu = 283.00296 * u.deg

Sun.frame = ICRS()
#sss=classical.ClassicalState(Sun, a, ecc, inc, raan, argp, nu,  epoch=t)
ss = Orbit.from_classical(Sun, a, ecc, inc, raan, argp, nu,  epoch=t)
ss._frame = ICRS()
print(ss.plane, ss.v)
