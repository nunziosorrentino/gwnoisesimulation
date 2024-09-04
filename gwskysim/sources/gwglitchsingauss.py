import numpy as np
from numpy.random import uniform

from gwskysim.sources.gwnoisesource import GWNoiseSource
from gwskysim.sources.gwpopulation import GWPopulation
from gwskysim.sources.gwsignal import GWSignal

class GWGlitchSinGauss(GWNoiseSource):

    def __init__(self,name,t0,tau,f0,SNR=None):
        GWNoiseSource.__init__(self, name, t0 ,tau,SNR)
        self.f0 = f0



    def __call__(self, times_array):
        times_array = times_array[times_array <= self.t0 + 5 * self.tau]
        times_array = times_array[times_array >= self.t0 - 5 * self.tau]
        if np.array(times_array).size==0:
            return None 

        Q = self.tau * ((2 ** 0.5) * np.pi * self.f0)


        if self.SNR is None:
        #if self.SNR is None and self.detector is None :
            #print("no SNR provided--> random amplitude")
            #print("SNR",self.SNR)
            #print("asd",self.detector)
            h0 =  uniform(10**-23,10**-21)
            aux =  np.exp(-((times_array - self.t0) ** 2) / (2 * self.tau ** 2))
            aux2= h0*np.sin(2 * np.pi * self.f0 * (times_array - self.t0)) * aux
            return aux2,times_array[0]

        else:
            #SNR,f_asd_tab,asd_tab
            #print("computing with the provided SNR")
            h0 = 10 ** (-21)

            aux = h0 * np.exp(-((times_array - self.t0) ** 2) / (2 * self.tau ** 2))
            aux2= h0*np.sin(2 * np.pi * self.f0 * (times_array - self.t0)) * aux

            snr_test=GWNoiseSource.SNR(aux2,times_array,self.f_asd_tab,self.asd_tab)
            #print("snr test sin gaus",snr_test)
            AUX = (self.SNR / snr_test) * aux2

            return AUX,times_array[0]

