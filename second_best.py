from scipy.optimize import fsolve
import scipy.optimize as opt
from math import atan, tan, asin, cos, sin, atan2,sqrt
import numpy as np
def f(variables) :
    (ro0,r) = variables
    p= 0.8055337998759792
    q=0.7479822145017828
    cc=0.8729375973678445
    r2=0.9849409084612417
#    ro0=0.6054236
    second_eq = -ro0+p-q/(r**3)
    first_eq = -(r*r)+(p-q/(r**3))**2+2*cc*(p-q/(r**3))+r2
    return [first_eq, second_eq]

solution = opt.fsolve(f, (0.1,10) )
r=solution[1]
ro0=solution[0]

print('ro=',ro0,'r=',r)
