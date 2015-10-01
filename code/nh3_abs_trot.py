#/usr/bin/env python

import numpy as np
import scipy.constants
import pylab as pl
import math
from numpy import arange,array,ones,linalg
#from pylab import plot,show
#from matlablib import plot,show




#Input constants
#IMPARA A PRENDERLE DA COMANDO
boltzk=1.3806503e-16 #cm2 g s-2 K-1
planckh=6.626068e-27# cm2 g/s
vc = 299792.458 #!in km/s
#
# ammonia rotational constants (Hz)
B=298117.0e6
C=186726.0e6
#B=9.44430 #in cm-1
#C=6.19600

# flux calibration correction factor
tentencorrection =  0.346 / 0.258


#Start providing the measured parameters
#
#0) Peak Flux; 1) Linewidth; 2) Integrated Flux; 3) tau

lines_e2w = np.array([[0, 0, 0, 0],  #oneone
                    [0, 0, 0, 0],  #twotwo
                    [0, 0, 0, 0],  #threethree
                    [0, 0, 0, 0],  #fourfour
                    [0, 0, 0, 0],  #fivefive
                    [-0.336,7.3,-2.61,120], # sixsix
                    [-0.310,5.9,-1.95,45], # sevenseven
                    [0, 0, 0, 0],  #eight
                    [-0.230,5.6,-1.38,18], # nine
                    [-0.150*tentencorrection, 4.4, -0.7*tentencorrection, 0], # ten
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [-0.072, 3.9, -0.3, 0], # thirteen
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0]])
notau_list_e2w = [10,13]
Fcont_e2w = [0, 0, 0, 0, 0, 0.362, 0.362, 0, 0.346, 0.258*tentencorrection, 0, 0, 0.369, 0, 0]

#e2e is thick!!!
# lines_e2e = np.array([[0, 0, 0, 0],  #oneone
#                     [0, 0, 0, 0],  #twotwo
#                     [0, 0, 0, 0],  #threethree
#                     [0, 0, 0, 0],  #fourfour
#                     [0, 0, 0, 0],  #fivefive
#                     [0.432, 10.6, 4.89, 39], # sixsix
#                     [0.3210, 8.6, 2.82, 14], # sevenseven
#                     [0, 0, 0, 0],  #eight
#                     [0.099, 9.0, 0.94, 0.55], # nine (upper limit on tau ~3 b/c line flattening becomes severe above this); tau=0.55 for tex=175
#                     [0.059*tentencorrection, 13, 0.82*tentencorrection, 0.1], # ten (0.1 is arbitrary)
#                     [0, 0, 0, 0],
#                     [0, 0, 0, 0],
#                     [0.069, 5.1, 0.38, 0.34], # thirteen (0.1 is arbitrary) (tau=0.34 for tex=175)
#                     [0, 0, 0, 0],
#                     [0, 0, 0, 0],
#                     [0, 0, 0, 0],
#                     [0, 0, 0, 0]])
# notau_list_e2e = []
# Fcont_e2e = [0]*14

lines_m = lines_e2w
notau_list = notau_list_e2w
Fcont = Fcont_e2w


#                    
Fpeak = abs(lines_m[:,0])
Dv = lines_m[:,1]
Fint = abs(lines_m[:,2])
tau = lines_m[:,3]
print Fpeak,Dv,Fint, tau
#
#
#------------------------------------------------------------------------------------
#Let's now calculate the temperatures/densities for transtions w/o opacities estimates.
#
#notau_list =[10,12,13,14] 
#notau_list =[6,7,9] 
for j1 in notau_list:
    j1 = int(j1-1)  #Account for the fact that the 66 line corresponds to index=5
    tau[j1] = -math.log( 1- ( Fpeak[j1] / Fcont[j1] ) )


#
#
# All metastable transitions. 0) frequencies in GHz, 1) J, 2) El (in K); 3) Snu^2 
lines_nh3 = np.array([[23.69450, 1,  24,     0],
                      [23.72263, 2,  65,     0],
                      [23.87013, 3, 124,     0],
                      [24.13942, 4, 201,     0],
                      [24.53299, 5, 296,     0],
                      [25.05596, 6, 408, 52.53],
                      [25.71514, 7, 538, 61.67],
                      [26.51896, 8, 686,     0],
                      [27.47794, 9, 852, 80.61],
                      [28.60475,10,1035, 90.00],
                      [29.91449,11,1236, 99.39],
                      [31.42494,12,1455,108.78],
                      [33.15684,13,1691,118.18],
                      [35.13428,14,1945,127.59],
                      [37.38513,15,2215,     0]])
freq = lines_nh3[:,0]
J = lines_nh3[:,1]
El = lines_nh3[:,2]
Snu2 = lines_nh3[:,3]
#
#
#
#------------------Here we calculate Trot from a Boltzmann Diagram--------------
#
#We use this formula that requires freq in GHz, Fpeak in Jy, angle in arcseconds
#        Nu / gu = 2.04 * Fpeak[val-1] * 1.0E20 / (theta_a * theta_b * Snu2[val-1] * freq[val-1]**3 )  


def t_d_tau(lines_list):
    """
    Estimates the rotational tempearature and column density 
    from a Boltzman Diagram and a least square fit. 
    It takes into account the opacity of the lines. 
    If there is not measurement of the opacity, uses the function tau to get one. 
    """
    Njk_over_Tex = []
    El_nh3 = []
    for val in lines_list:
        print tau[val-1]
        b = 1.61 * 1.0E14 * (J[val-1]+1)/(J[val-1]*(2*J[val-1]+1)) * tau[val-1] * Dv[val-1] / freq[val-1]
        if (J[val-1]==3 or J[val-1]==6  or J[val-1]==9 or J[val-1]==12 or J[val-1]==15):
            b=b/2.0
        Njk_over_Tex.append( b )
        El_nh3.append(El[val-1])
    y = np.log(Njk_over_Tex)
    A = np.array([El_nh3, np.ones(len(El_nh3))])
    m, c = np.linalg.lstsq(A.T,y)[0]
    Trot=(-1.0/m)
    print "Trot = {0}".format(Trot)
    if np.isnan(Trot):
        import ipdb; ipdb.set_trace()

    #calculating partition function Q_rot; use egn 15.48 tools of radio astronomy
    # = sqrt( (pi (kT)^3) / h^3 B^2 C )
    Q_rot= math.sqrt( math.pi * (boltzk*Trot)**3 / ( planckh**3 * B**2 * C ) )
    print "Q = {0}".format(Q_rot)
    # read off density from y axis, it is Result[0]=ln(N_T/Q_rot) 
    N_over_Tex = Q_rot * math.exp( c )

    p, = pl.semilogy(El_nh3, np.exp(y), 's', label='\n$T={0:0.1f}$\n$N(\\mathrm{{NH}}_3)={1:0.2g}$'.format(Trot, N_over_Tex))
    pl.plot(np.linspace(0,1800), np.exp(m*np.linspace(0,1800) + c), color=p.get_color())
    pl.xlabel('$E_u$ [K]')
    pl.ylabel("$N_u / g$ [cm$^{-2}$]")
    pl.legend(loc='best')


    return  Trot, N_over_Tex
#
#
#lines_list = [6,7]
#T_679, N_over_Tex_679 = t_d_tau(lines_list)
#print "For 66-77-99  the Rotational temperature is: ", T_679, "and the column density is:", N_over_Tex_679
#
#lines_list = [6,7,9]
#T_679, N_over_Tex_679 = t_d_tau(lines_list)
#print "For 66-77-99  the Rotational temperature is: ", T_679, "and the column density is:", N_over_Tex_679
##
#lines_list = [6,9]
#T_69, N_over_Tex_69 = t_d_tau(lines_list)
#print "For 66-99  the Rotational temperature is: ", T_69, "and the column density is:", N_over_Tex_69
#
#
#Let's now calculate the temperatures/densities for transtions w/o opacities estimates.
#
#T_10121314, rho_10121314 = t_d_tau(notau_list)
#print "For 1010-1212-1313-1414 the Rotational temperature is: ", T_10121314, "and the column density is:", rho_10121314
#a = [10,13,14] 
#T_101314, rho_101314 = t_d_tau(a)
#print "For 1010-1313-1414    the Rotational temperature is: ", T_101314, "and the column density is:", rho_101314
#a = [13,14]
#T_1314, rho_1314 = t_d_tau(a)
#print "For 1313-1414         the Rotational temperature is: ", T_1314, "and the column density is:", rho_1314

pl.figure(2).clf()
T,N = t_d_tau([6,7,9,13,])
from astropy import units as u

# 2000 au comes from an integration region 0.8 by 0.9"
radius = 2000*u.au # was: 0.0078*u.pc
radius = 1377*u.au # from gaussian fitting - deconvolved source size

Xnh3 = 1e-7
col = ((N*u.cm**-2)/Xnh3 / (radius)).to(u.cm**-3)
print col, np.log10(col.value)
mass = ((N*u.cm**-2)/Xnh3 * (2*np.pi*(radius)**2) * (2.8*u.Da)).to(u.M_sun)
print mass, np.log10(mass.value)

pl.savefig('figures/absorption_physical_parameter_fit_hardcodedtau.png')


pl.draw()
pl.show()
