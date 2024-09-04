# classe molto generale che eredita i metodi di GWSource
import os
import numpy as np
from gwskysim.sources.gwsource import GWSource
from gwskysim.sources.gwsignal import GWSignal
from pycbc.types.timeseries import TimeSeries
from pycbc.detector import Detector
from gwskysim.sources.gwsignal import GWSignal
#from gwskysim.sources.gwpopulation import GWPopulation
from gwskysim import GWSKYSIM_CONF 

asd_dict = {"V1": os.path.join(GWSKYSIM_CONF, "V1_O3.txt"),
            "H1": os.path.join(GWSKYSIM_CONF, "H1_O3.txt"),
            "L1": os.path.join(GWSKYSIM_CONF, "L1_O3.txt")
            }

class GWNoiseSource(GWSource):

    def __init__(self, name, t0, tau, SNR=None):
        GWSource.__init__(self, name)
        self.t0 = t0
        self.tau = tau
        self.SNR = SNR
        global snr
        snr=SNR
        #self.det=None
        #print("SNR",SNR)

    def set_dect(self,dct):

        if dct in ["V1", "L1", "H1"]:
            f_asd_tabV, asd_tabV = np.loadtxt(asd_dict[dct], unpack=True)           
            self.SNR=snr
            self.f_asd_tab = f_asd_tabV
            self.asd_tab = asd_tabV
            #print("get dect2", dct)

        else:
            self.SNR = None
            self.f_asd_tab = None
            self.asd_tab = None
            #print("no dect")

    @staticmethod
    def SNR(signal_amplitude,signal_times,f_asd_tab,asd_tab):

        samples = len(signal_times)
        dt = np.abs(signal_times[1] - signal_times[0])
        #print("dt", dt)

        f = rfftfreq(samples, dt)
        df = np.abs(f[2] - f[1])
        #print("f:",f)

        # print("frequency",f)

        fft_signal = (np.fft.rfft(signal_amplitude)) / samples
        #print(np.abs(fft_signal))
        f_max_asd = np.amax(f_asd_tab)
        f_min_asd = np.amin(f_asd_tab)

        F_max_sig = np.amax(f)
        F_min_sig = np.amin(f)
        #print("---------------",F_max_sig)
        asd_interp = interp1d(f_asd_tab, asd_tab)

        if F_max_sig < f_min_asd:
            #print("ERROR: signal not in the detector frequency band")
            #print("F_max_signal", F_max_sig)
            #print("f_min_asd", f_min_asd)

            SNR=0

        if F_min_sig > f_max_asd:

            #print("ERROR: signal not in the detector frequency band")
            #print("F_min_signal", F_min_sig)
            #print("f_max_asd", f_max_asd)
            SNR=0


        else:
            """
            print("F_max_signal", F_max_sig)
            print("F_min_signal", F_min_sig)
            print("f_min_asd", f_min_asd)
            print("f_max_asd", f_max_asd)
            """


            aux = np.where(f >= f_min_asd, f, 0 * f)
            aux2 = np.nonzero(aux)
            index = aux2[0][0]  # indice da cui paartire a campionare il segnale
            #print("indexxxxxxxxxxxx",index)
            if F_max_sig < f_max_asd:

                f_aux = np.linspace(f_min_asd, F_max_sig,  len(f[index:]))
                #print("len signal", len(f[index:]))
                #print("len2",len(f_aux))

                integrand = 4 * (np.absolute(fft_signal[index:]) ** 2) / (np.absolute(asd_interp(f_aux))) ** 2
                SNR = integrate.simps(integrand, dx=df)
                #plt.plot((np.absolute(fft_signal[index:]) ** 2))
                #plt.xscale("log")
                #plt.yscale("log")
                #print("---------------------------",SNR**0.5)

            else:

                aux3 = np.where(f < f_max_asd, 0 * f, f)
                aux4 = np.nonzero(aux3)
                index2 = aux4[0][0]

                f_aux = np.linspace(f_min_asd, f_max_asd, len(f[index:index2]))
                #print("len signal", len(f[index:index2]))

                integrand = 4 * (np.absolute(fft_signal[index:index2]) ** 2) / (np.absolute(asd_interp(f_aux))) ** 2
                SNR = integrate.simps(integrand, dx=df)

        return SNR ** 0.5 


    def _import_template(self, times_array, dect=None):
        
        self.set_dect(dect)
        
        if len(times_array) < 2:
            raise ValueError('Time Series must contain at least two sample.')
        if self(times_array) is None:
            return self(times_array) 
        try:
            h_s,time0 = self(times_array)
        except ValueError:
            h_s=self(times_array)
            time0=times_array[0]
        m_sample_rate = round(1/(times_array[1]-times_array[0]))
        dt_=len(h_s)/m_sample_rate
        h=GWSignal(self.name,time0,dt_,m_sample_rate,h_s,detector=dect)

        return h

    def generate_signal(self, m_start_time, m_duration, m_sample_rate,
                        detectors=['H1', 'L1', 'V1'], approximant=None, 
                        lower_freq=40):

        signal_projected_all_det=[]

        for ds in detectors:
            #print("i",ds)
            self.det=ds
            times_array = np.linspace(m_start_time, m_start_time + m_duration, 
                                      np.int(m_duration * m_sample_rate))
            h = self._import_template(times_array, dect=self.det)
            if h is None:
                return h

            signal_one_det = GWSignal(self.name,
                                    start_time=h.start_time,
                                    duration=h.duration,
                                    sample_rate=m_sample_rate,
                                    sampled_strain=h.sampled_strain,
                                    detector=ds)

            signal_projected_all_det.append(signal_one_det)
        #print("len!",len(signal_projected_all_det))
        return signal_projected_all_det

       # print((h))

        #return h


