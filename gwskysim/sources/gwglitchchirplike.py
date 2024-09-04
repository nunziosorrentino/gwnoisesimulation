import numpy as np

from gwskysim.sources.gwnoisesource import GWNoiseSource


class GWGlitchChirpLike(GWNoiseSource):        #        m1  m2 solar mass

    def __init__(self, name, t0, m1,m2,SNR=None):  #tau is redundant and it is set to 0
        GWNoiseSource.__init__(self, name, t0, 0,SNR)
        self.m1=m1*1.989*10**(30)
        self.m2=m2*1.989*10**(30)


    def __call__(self, times_array):
        times_array = times_array[times_array <= self.t0 ]
        times_array = times_array[times_array >= self.t0 - 5]
        if np.array(times_array).size==0:
            return None 

        if self.SNR is None:
        #if self.SNR is None and self.detector is None :

            #print("random amplitude")
            h0 = np.random.uniform(10 ** (-22), 5 * 10 ** (-21))

            chirp_mass=((self.m1*self.m2)**(3./5.))/((self.m1+self.m2)**(1./5))

            l = len(times_array)

            AuX2 = np.zeros(l)
            tau=np.zeros(l)
            phase=np.zeros(l)

            for i in range(l):
                # print(times_array[i])

                if times_array[i] >self.t0-0.2 and times_array[i] < self.t0+0.0009:

                    tau[i]=np.abs(self.t0-times_array[i])
                    phase[i]=(-2*((5*6.67*(10**(-11))*chirp_mass)/(3*10**8)**3)**(-5./8))*(tau[i])**(5/8)



                    AuX2[i] = h0*np.cos(phase[i])*(1/tau[i])**(1/4)

            return AuX2,times_array[0]

        else :
            #print("computing SNR whit provided sensitivity curves ")
            ### SNR, f_asd_tab, asd_tab
            h0 = 10 ** (-20)

            chirp_mass = ((self.m1 * self.m2) ** (3. / 5.)) / ((self.m1 + self.m2) ** (1. / 5))

            l = len(times_array)

            AuX2 = np.zeros(l)
            tau = np.zeros(l)
            phase = np.zeros(l)

            for i in range(l):
                # print(times_array[i])
                # di default il chirp dura 5 secondi, non è molto realistico..., meglio implementare formula 4.21 del Maggiore vol1 .

                if times_array[i] >self.t0-5 and times_array[i] < self.t0:
                    tau[i] = np.abs(self.t0 - times_array[i])
                    phase[i] = (-2 * ((5 * 6.67 * (10 ** (-11)) * chirp_mass) / (3 * 10 ** 8) ** 3) ** (-5. / 8)) * (
                    tau[i]) ** (5 / 8)

                    AuX2[i] = h0 * np.cos(phase[i]) * (1 / tau[i]) ** (1 / 4)
                """
                if times_array[i] > self.t0  and times_array[i] < self.t0+0.001:
                    tau[i] = np.abs(self.t0 - times_array[i])
                    phase[i] =(-2 * ((5 * 6.67 * (10 ** (-11)) * chirp_mass) / (3 * 10 ** 8) ** 3) ** (-5. / 8)) * (
                    tau[i]) ** (5 / 8)


                    AuX2[i] = 1*h0 * np.cos(phase[i]) * (1 / tau[i]) ** (1 / 4)
                """

            snr_test = GWNoiseSource.SNR(AuX2, times_array, self.f_asd_tab,self.asd_tab)
            #print("aaaaaaa",snr_test)
            return (self.SNR / snr_test) * AuX2,times_array[0]
            #return  AuX2, times_array[0]


