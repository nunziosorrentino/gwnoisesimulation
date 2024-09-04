import pandas as pd
import h5py
import numpy as np

from astropy import units as u
from astropy.constants import c, G, M_sun

def pc_to_Msun(dist):
    """
    """
    #dist = dist*((27/14)*(10**16)) 
    dist = dist * u.parsec 
    dist = dist.to(u.m).value
    dist = dist*(c.value**2)/(G.value*M_sun.value)
    #print(c, G, M_sun)
    return dist
    
def Msun_to_sec(mass):
    """
    """
    t_ = mass*(G.value*M_sun.value)/(c.value**3)
    #print(c, G, M_sun)
    return t_   

def random_from_distribution(distribution, val_min, val_max, length):
    """
    Generate a list of random numbers from a given distribution in a given interval
    :param distribution: at the moment must be 'uniform','cos' or 'sin' (means uniform in cos/sin: use 'sin' for
    declination if you want a uniform sky distribution)
    :param val_min: minimum value
    :param val_max: maximum value
    :param length: how many random numbers you want in the list
    :return: a list of random numbers
    """
    if distribution == 'uniform':
        random_sample = np.random.uniform(low=val_min, high=val_max, size=(length,))
        return random_sample
    if distribution == 'sin':
        random_sample = []
        i = 0
        while i < length:
            random_num = np.arcsin(2 * np.random.uniform() - 1)
            random_num = np.rad2deg(random_num)
            if val_min < random_num < val_max:
                i = i + 1
                random_sample.append(random_num)
        return np.array(random_sample)
    if distribution == 'cos':
        random_sample = []
        i = 0
        while i < length:
            random_num = np.arccos(2 * np.random.uniform() - 1)
            random_num = np.rad2deg(random_num)
            if val_min < random_num < val_max:
                i = i + 1
                random_sample.append(random_num)
        return np.array(random_sample)

    if distribution == "delta1":
       #equally spaced between val_max e val_min
        aux = val_max - val_min
        aux2 = aux / length
        vector = np.zeros(length)
        #print("len",len(vector))
        for i in range(length):
            #print("i", i)
            if i == 0:
                vector[0] = val_min + aux2 / 2
            if i == length :
                vector[length ] = val_max - aux2 / 2
            else:
                vector[i] = vector[0] + i * aux2
        return vector

    if distribution == "growing_random":

        aux = val_max - val_min
        aux2 = aux / length
        vector = np.zeros(length)
        #print("len",len(vector))
        for i in range(length):
            #print("i", i)
            if i == 0:
                vector[0] = np.random.uniform(low=val_min,high=val_min + aux2 / 2)
            elif i == length :
                vector[length ] = np.random.uniform(low=val_max - aux2 / 2,high=val_max)
            else:
                centered_f=val_min+i*aux2
                vector[i] = np.random.uniform(low=centered_f-aux2/2 ,high=centered_f+aux2/2)
        return vector
        
    if distribution == "power3":
        val_dist = val_max - val_min
        random_sample = np.random.power(4, size=length)
        return random_sample*val_dist + val_min    

