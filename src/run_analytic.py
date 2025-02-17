import numpy as np
import sys
from analytic import *

print('Example usage: python main.py theta')

# setting up model parameters
a1 = -0.5 # bottom depth of the anisotropic layer
a2 = -0.1 # top depth
mid_layer_depth = (a1 + a2)/2
w  = 1
ux0 = 1
theta = float(sys.argv[1]) # 
n1 = np.cos((theta+90)/180*np.pi)
n2 = np.sin((theta+90)/180*np.pi)
es = 0.1 # weak viscosity
e  = 1 # strong viscosity
dd = 0.1 # spatial grids along the depth profile.
mid_layer_depth = (a1 + a2)/2
id = round((mid_layer_depth-(-w))/dd+1) # counting from -w to 0, idx 8 is inside a1 ~ a2 depths of the anisotropic zone.

# compute
d, sig11, sig12, sig22, str11, str12, str22, u1, p = analytic(a1, a2, w, ux0, n1, n2, es, e, dd)
smax, smin, nx0, ny0, nx1, ny1, J2 = calc_principal(sig11[id-1], sig22[id-1], sig12[id-1])
srmax, srmin, srnx0, srny0, srnx1, srny1, srJ2 = calc_principal(str11[id-1], str22[id-1], str12[id-1])
srmax_iso, srmin_iso, srnx0_iso, srny0_iso, srnx1_iso, srny1_iso, srJ2_iso = calc_principal(str11[2], str22[2], str12[2])
theta_sigma1 = np.arctan(ny0/nx0)/np.pi*180
theta_eps1 = np.arctan(srny0/srnx0)/np.pi*180

#print(smax, smin, nx0, ny0, nx1, ny1, J2)
#print(srmax, srmin, srnx0, srny0, srnx1, srny1, srJ2)
print(f"Normal director to weak anisotropy direction rotates from y+ at {theta:.2f} degs")
print(f"The angle of sigma1 to x+ is {theta_sigma1:.2f}")
print(f"The angle of eps1 to x+ is {theta_eps1:.2f}")
print(f"Strain localization is {srJ2/srJ2_iso:.2f}")
print(f"Pressure inside ani and iso layers are {p[id-1]:.2f} and {p[1]:.2f}")
