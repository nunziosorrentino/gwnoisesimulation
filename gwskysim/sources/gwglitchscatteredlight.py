import numpy as np

from gwskysim.sources.gwnoisesource import GWNoiseSource

class GWGlitchScatteredLight(GWNoiseSource): #scattered light glitch

    def __init__(self,name,t0,tau, f0,SNR=None):
        GWNoiseSource.__init__(self, name, t0 ,tau,SNR)
        self.f0 = f0
        self.multipl=12

    def __call__(self, times_array):
        times_array = times_array[times_array <= self.t0 + 8 * self.tau]
        times_array = times_array[times_array >= self.t0 - 8 * self.tau]
        if np.array(times_array).size==0:
            return None         
        if self.SNR is None:
        # if self.SNR is None and self.detector is None :
            #print("no SNR --->random amplitude")
            h0 = np.random.uniform(10 ** (-22),  5* 10 ** (-21))
            phi =2.* np.pi * self.f0*np.abs(times_array-self.t0)*(1.-0.5*(times_array-self.t0)**2.)
            return h0*np.sin(phi)*np.exp(-(times_array-self.t0)**2/(2*self.tau)),times_array[0]

        else:
            #print("computing with the given snr ",self.SNR)
            h0 = 10 ** (-10)

            phi = 2. * np.pi * self.f0 * np.abs(times_array - self.t0) * (1. - 0.5 * (times_array - self.t0) ** 2.)
            aux = h0 * np.sin(phi) * np.exp(-(times_array - self.t0) ** 2 / (2 * self.tau))
            snr_test = GWNoiseSource.SNR(aux, times_array, self.f_asd_tab, self.asd_tab)
            # print("self.snr",aux)
            # print("snr test",snr_test)
            return (self.SNR / snr_test) * aux,times_array[0]



            """
            leng=self.tau/(times_array[1]-times_array[0])
            times_array_reduced=times_array[0:np.int(12.*leng)] #finestra
            t0=times_array[np.int(6.*leng)]

            phi = 2. * np.pi * self.f0 * np.abs(times_array_reduced - t0) * (1. - 0.5 * (times_array_reduced - t0) ** 2.)
            aux = h0 * np.sin(phi) * np.exp(-(times_array_reduced - t0) ** 2 / (2 * self.tau))

            snr_test = GWNoiseSource.SNR(aux, times_array_reduced, self.f_asd_tab, self.asd_tab)

            AUX= (self.SNR / snr_test) * aux
            AUX2= GWNoiseSource.gluer(AUX,times_array_reduced,self.t0,times_array)

            return AUX2
            """


"""
h0 = 10 ** -21
phi = 2. * np.pi * self.f0 * np.abs(times_array - self.t0) * (1. - 0.5 * (times_array - self.t0) ** 2.)
aux = h0 * np.sin(phi) * np.exp(-(times_array - self.t0) ** 2 / (2 * self.tau))

snr_test = GWNoiseSource.SNR(aux, times_array, self.f_asd_tab, self.asd_tab)
# print("self.snr",aux)
# print("snr test",snr_test)
return (self.SNR / snr_test) * aux

"""
