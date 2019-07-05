
# these functions replicate figures in van Genuchten, M.Th., 1981. Analytical solutions for chemical transport with simultaneous adsorptions, zero-order production and first-order decay. Journal of Hydrology 49. 213--233.

import chemtran

D = 0.00375 # chemical dispersion [0.1,5.0] m²/day
R = 3.0 # retardation factor [1.0,25.0]
mu = 0.25 
v = 0.25 # interstitial/pore-water velocity [m/day] (=q/theta); theta=volumetric moisture content ~ porosity at saturated
# k = 0.1 # first-order decay coefficient [-7.0,-0.02]; zero when not needed (alpha in van Genuchten, 1981)

ct = chemtran.chemtran(D,v,R,mu)

Ci = 0.0 # initial concentration [meq/L]
C0 = 1.0 # influx concentration [meq/L]
t0 = 5.0 # i.e., a five-day pulse


# # Case A1 (steady pulse)
# # Test1: fig.2
# gs = [0.75, 0.5, 0.25, 0.0]
# for x in range(0, 101):
#     xx = x/100
#     ln = str(xx)
#     for g in gs: # gamma: zero-order liquid phase source term meq/L/day
#         ln += "," + str(ct.cxtA1(xx, 7.5, C0, t0, g))
#     print(ln)

# # Test2: fig.3
# mus = [0.0, 0.25, 0.5, 1.0]
# for x in range(0, 101):
#     xx = x/100
#     ln = str(xx)
#     for mu in mus:
#         ct.mu = mu
#         ln += "," + str(ct.cxtA1(xx, 7.5, C0, t0, 0.5))
#     print(ln)


# # Case A2 (steady-state)
# ct.mu = 0.5
# Cb = 0.0
# # Test3: fig.4 (part1)
# for x in range(0, 101):
#     xx = x/100
#     ln = str(xx)
#     ln += "," + str(ct.cxA2(xx, C0, 0.25))
#     print(ln)

# # Test3: fig.4 (part2)
# for x in range(0, 101):
#     xx = x/100
#     ln = str(xx)
#     ln += "," + str(ct.cxiA3(xx, Cb, 0.25))
#     print(ln)    

# # Test3: fig.4 (part3)
# ts = [2.5,5.0,7.5,10.0]
# for x in range(0, 101):
#     xx = x/100
#     ln = str(xx)
#     for t in ts:
#         ln += "," + str(ct.cxtA3(xx, t, Cb, C0, t0, 0.25))
#     print(ln) 