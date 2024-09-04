from gwskysim.sources.gwnoisesource import GWNoiseSource

import numpy as np
#from spectrum.arma import arma_estimate, arma2psd
#from statsmodels.tsa.arima_process import arma_generate_sample

from scipy import signal

#funzioni che genero rumore bianco gaussiano e non (ARMA) con lunghezza pari a times_array dei glitch


class White_gaussian_noise(GWNoiseSource): #white gaussian noise with zero mean and std equals to 1

    def __init__(self,name):
        GWNoiseSource.__init__(self, name, t0=0,tau=0)


    def __call__(self, times_array,dect=None):

        l=len(times_array)
        h0 = np.random.uniform(1 * 10 ** (-21), 5 * 10 ** (-21))
        s = np.random.normal(0, 1, l)
        return s*h0


#class White_noise_ARMA(GWNoiseSource): #white noise

#    def __init__(self,name):
#        GWNoiseSource.__init__(self, name, t0=0,tau=0)
#
#
#    def __call__(self,times_array,dect=None):
#
#        N=len(times_array)
#        arparams = np.array([.75, -.25])
#        maparams = np.array([.65, .35])

#        arparams = np.r_[1, -arparams]
#        maparams = np.r_[1, maparams]
#        nobs = N
#        y = arma_generate_sample(arparams, maparams, nobs)

#        h0 = np.random.uniform(1 * 10 ** (-21), 5 * 10 ** (-21))

#        return h0*(y / np.amax(y))

