# ATTENTION
# This is a script that group a set of time domain waveform functions for the close-encounters modelling 
# see N. Loutrel 2019 for details

import numpy as np
import pickle
import os
from gwskysim import GWSKYSIM_WD

cos_file = open(os.path.join(GWSKYSIM_WD, "models", "cosarray.pickle"), "rb")
cos_array = pickle.load(cos_file)
cos_file.close()
sin_file = open(os.path.join(GWSKYSIM_WD, "models", "sinarray.pickle"), "rb")
sin_array = pickle.load(sin_file)
sin_file.close()

# Extract function:

def sintable(angle):
    index = int(199*angle/(2*np.pi))
    return sin_array[index] 
    
def costable(angle):
    index = int(199*angle/(2*np.pi))
    return cos_array[index]     

# Easy functions

def th_ci_si(incl):
    return 3 + costable(incl)**2 - sintable(incl)**2

def cb_sb(pol):
    return costable(pol)**2 - sintable(pol)**2

############################# 
#                           #
# cosin-like plus functions #
#                           #
#############################

def cp_0_0(ph, incl, pol):
    first_ = 8*ph/(1 + np.power(ph, 2))
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl))
    
def cp_0_1(ph, incl, pol):
    first_ = 9*costable(pol)**2*th_ci_si(incl)
    second_ = -9*(1 + costable(incl)**2*(-1 + sintable(pol)**2) + sintable(incl)**2 - sintable(pol)**2*(-3 + sintable(incl)**2))
    third_ = (16*costable(pol)*sintable(pol)*th_ci_si(incl)*ph)/(1 + np.power(ph, 2))
    return np.array((first_ + second_ + third_)/20)

def cp_0_2(ph, incl, pol):
    first_ = -9*costable(incl)**2*(3 + 265*sintable(pol)**2)
    second_ = 2385*costable(pol)**2*th_ci_si(incl)
    third_ = 2385*sintable(pol)**2*(-3 + sintable(incl)**2)
    fourth_ = 27*(1 + sintable(incl)**2) 
    fifth_ = (3392*costable(pol)*sintable(pol)*th_ci_si(incl)*ph)/(1 + np.power(ph, 2))
    return np.array((first_ + second_ + third_ + fourth_ + fifth_)/2800)

def cp_1_0(ph, incl, pol):
    first_ = 1/np.power(1 + np.power(ph, 2), 1.5)
    second_ = -15*sintable(pol)**2 + sintable(incl)**2 + 5*sintable(pol)**2*sintable(incl)**2 + np.power(ph, 2) 
    third_ = 3*sintable(pol)**2*np.power(ph, 2) + sintable(incl)**2*np.power(ph, 2) - sintable(pol)**2*sintable(incl)**2*np.power(ph, 2)  
    fourth_ = -costable(pol)**2*th_ci_si(incl)*(-5 + np.power(ph, 2))
    fifth_ = costable(incl)**2*(-1 - np.power(ph, 2) + sintable(pol)**2*(-5 + np.power(ph, 2)))
    return first_*(1 + second_ + third_ + fourth_ + fifth_)

def cp_1_1(ph, incl, pol):
    first_ = -0.3/np.power(1 + np.power(ph, 2), 1.5)
    second_ = -15*sintable(pol)**2 + sintable(incl)**2 + 5*sintable(pol)**2*sintable(incl)**2 + np.power(ph, 2) 
    third_ = -21*sintable(pol)**2*np.power(ph, 2) + sintable(incl)**2*np.power(ph, 2) + 7*sintable(pol)**2*sintable(incl)**2*np.power(ph, 2)  
    fourth_ = costable(pol)**2*th_ci_si(incl)*(5 + 7*np.power(ph, 2))
    fifth_ = -costable(incl)**2*(1 + np.power(ph, 2) + sintable(pol)**2*(5 + 7*np.power(ph, 2)))
    return np.array(first_*(1 + second_ + third_ + fourth_ + fifth_))

def cp_1_2(ph, incl, pol):
    first_ = 7.1429e-4/np.power(1 + np.power(ph, 2), 1.5)
    second_ = 3585*sintable(pol)**2 - 111*sintable(incl)**2 - 1195*sintable(pol)**2*sintable(incl)**2 - 111*np.power(ph, 2) 
    third_ = 6915**sintable(pol)**2*np.power(ph, 2) - 111*sintable(incl)**2*np.power(ph, 2) - 2305*sintable(pol)**2*sintable(incl)**2*np.power(ph, 2)  
    fourth_ = -5*costable(pol)**2*th_ci_si(incl)*(239 + 461*np.power(ph, 2))
    fifth_ = costable(incl)**2*(111*(1 + np.power(ph, 2)) + 5*sintable(pol)**2*(239 + 461*np.power(ph, 2)))
    return np.array(first_*(-111 + second_ + third_ + fourth_ + fifth_))

def cp_2_0(ph, incl, pol):
    first_ = -2/(1 + np.power(ph, 2))
    return np.array(first_*cb_sb(pol)*th_ci_si(incl))

def cp_2_1(ph, incl, pol):
    first_ = -0.2/(1 + np.power(ph, 2))
    return np.array(first_*cb_sb(pol)*th_ci_si(incl))

def cp_2_2(ph, incl, pol):
    first_ = 0.00714/(1 + np.power(ph, 2))
    second_ = 96*sintable(pol)**2 + sintable(incl)**2 - 32*sintable(pol)**2*sintable(incl)**2 + np.power(ph, 2) 
    third_ = 3**sintable(pol)**2*np.power(ph, 2) + sintable(incl)**2*np.power(ph, 2) - sintable(pol)**2*sintable(incl)**2*np.power(ph, 2)  
    fourth_ = -costable(pol)**2*th_ci_si(incl)*(32 + np.power(ph, 2))
    fifth_ = costable(incl)**2*(-1 - np.power(ph, 2) + sintable(pol)**2*(32 + np.power(ph, 2)))
    return np.array(first_*(1 + second_ + third_ + fourth_ + fifth_))

def cp_4_2(ph, incl, pol):
    first_ = 0.0429/(1 + np.power(ph, 2))
    return np.array(first_*cb_sb(pol)*th_ci_si(incl))

########################### 
#                         #
# sin-like plus functions #
#                         #
###########################

def sp_1_0(ph, incl, pol):
    first_ = 8/(1 + np.power(ph, 2))
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl))

def sp_1_1(ph, incl, pol):
    first_ = -8/(5*(1 + np.power(ph, 2)))
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl))

def sp_1_2(ph, incl, pol):
    first_ = 4/(175*(1 + np.power(ph, 2)))
    second_ = -179 + 26*(1 + np.power(ph, 2))
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl)*second_)

def sp_2_0(ph, incl, pol):
    first_ = 12/np.power(1 + np.power(ph, 2), 1.5)
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl))

def sp_2_1(ph, incl, pol):
    first_ = -2/(5*np.power(1 + np.power(ph, 2), 1.5))
    second_ = 19 + 22*np.power(ph, 2)
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl)*second_)

def sp_2_2(ph, incl, pol):
    first_ = -1/(350*np.power(1 + np.power(ph, 2), 1.5))
    second_ = 2251 + 2026*np.power(ph, 2)
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl)*second_)

def sp_4_2(ph, incl, pol):
    first_ = 78/(35*np.power(1 + np.power(ph, 2), 1.5))
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl))

def sp_5_1(ph, incl, pol):
    first_ = 12/(5*(1 + np.power(ph, 2)))
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl))

def sp_5_2(ph, incl, pol):
    first_ = 36/(25*(1 + np.power(ph, 2)))
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl))

def sp_6_0(ph, incl, pol):
    first_ = -4/np.power(1 + np.power(ph, 2), 1.5)
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl))

def sp_6_1(ph, incl, pol):
    first_ = -2/(5*np.power(1 + np.power(ph, 2), 1.5))
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl))

def sp_6_2(ph, incl, pol):
    first_ = -37/(70*np.power(1 + np.power(ph, 2), 1.5))
    return np.array(first_*costable(pol)*sintable(pol)*th_ci_si(incl))

############################## 
#                            #
# cosin-like cross functions #
#                            #
##############################

def cc_0_0(ph, incl, pol):
    first_ = 16*ph/(1 + np.power(ph, 2))
    return np.array(first_*costable(incl)*cb_sb(pol))

def cc_0_1(ph, incl, pol):
    first_ = -2/(5*(1 + np.power(ph, 2)))
    second_ = -4*costable(pol)**2*ph + 4*sintable(pol)**2*ph + 9*costable(pol)*sintable(pol)*(1 + np.power(ph, 2))
    return np.array(first_*costable(incl)*second_)

def cc_0_2(ph, incl, pol):
    first_ = -53/(350*(1 + np.power(ph, 2)))
    second_ = -16*costable(pol)**2*ph + 16*sintable(pol)**2*ph + 45*costable(pol)*sintable(pol)*(1 + np.power(ph, 2))
    return np.array(first_*costable(incl)*second_)

def cc_1_0(ph, incl, pol):
    first_ = 1/(35*np.power(1 + np.power(ph, 2), 1.5))
    second_ = -1400 + 280*np.power(ph, 2)
    return np.array(first_*costable(pol)*sintable(pol)*costable(incl)*second_)

def cc_1_1(ph, incl, pol):
    first_ = 1/(35*np.power(1 + np.power(ph, 2), 1.5))
    second_ = 420 + 588*np.power(ph, 2)
    return np.array(first_*costable(pol)*sintable(pol)*costable(incl)*second_)

def cc_1_2(ph, incl, pol):
    first_ = 1/(35*np.power(1 + np.power(ph, 2), 1.5))
    second_ = 239 + 461*np.power(ph, 2)
    return np.array(first_*costable(pol)*sintable(pol)*costable(incl)*second_)

def cc_2_0(ph, incl, pol):
    first_ = 16/(1 + np.power(ph, 2))
    return np.array(first_*costable(pol)*sintable(pol)*costable(incl))

def cc_2_1(ph, incl, pol):
    first_ = 8/(5*(1 + np.power(ph, 2)))
    return np.array(first_*costable(pol)*sintable(pol)*costable(incl))

def cc_2_2(ph, incl, pol):
    first_ = 2/(35*(1 + np.power(ph, 2)))
    second_ = 32 + np.power(ph, 2)
    return np.array(first_*costable(pol)*sintable(pol)*costable(incl)*second_)

def cc_4_2(ph, incl, pol):
    first_ = -12/(35*(1 + np.power(ph, 2)))
    return np.array(first_*costable(pol)*sintable(pol)*costable(incl))

############################ 
#                          #
# sin-like cross functions #
#                          #
############################

def sc_1_0(ph, incl, pol):
    first_ = 16/(1 + np.power(ph, 2))
    return np.array(first_*costable(incl)*cb_sb(pol))

def sc_1_1(ph, incl, pol):
    first_ = -16/(5*(1 + np.power(ph, 2)))
    return np.array(first_*costable(incl)*cb_sb(pol))

def sc_1_2(ph, incl, pol):
    first_ = 8/(175*(1 + np.power(ph, 2)))
    second_ = -179 + 26*np.power(ph, 2)
    return np.array(first_*costable(incl)*cb_sb(pol)*second_)

def sc_2_0(ph, incl, pol):
    first_ = 24/np.power(1 + np.power(ph, 2), 1.5)
    return np.array(first_*costable(incl)*cb_sb(pol))

def sc_2_1(ph, incl, pol):
    first_ = -4/(5*np.power(1 + np.power(ph, 2), 1.5))
    second_ = 19 + 22*np.power(ph, 2)
    return np.array(first_*costable(incl)*cb_sb(pol)*second_)

def sc_2_2(ph, incl, pol):
    first_ = -1/(175*np.power(1 + np.power(ph, 2), 1.5))
    second_ = 2251 + 2026*np.power(ph, 2)
    return np.array(first_*costable(incl)*cb_sb(pol)*second_)

def sc_4_2(ph, incl, pol):
    first_ = 156/(35*np.power(1 + np.power(ph, 2), 1.5))
    return np.array(first_*costable(incl)*cb_sb(pol))

def sc_5_1(ph, incl, pol):
    first_ = 24/(5*(1 + np.power(ph, 2)))
    return np.array(first_*costable(incl)*cb_sb(pol))

def sc_5_2(ph, incl, pol):
    first_ = 72/(25*(1 + np.power(ph, 2)))
    return np.array(first_*costable(incl)*cb_sb(pol))

def sc_6_0(ph, incl, pol):
    first_ = -8/np.power(1 + np.power(ph, 2), 1.5)
    return np.array(first_*costable(incl)*cb_sb(pol))

def sc_6_1(ph, incl, pol):
    first_ = -4/(5*np.power(1 + np.power(ph, 2), 1.5))
    return np.array(first_*costable(incl)*cb_sb(pol))

def sc_6_2(ph, incl, pol):
    first_ = -37/(35*np.power(1 + np.power(ph, 2), 1.5))
    return np.array(first_*costable(incl)*cb_sb(pol))


