from gwskysim.sources.gwnoisesource import GWNoiseSource


import numpy as np
from numpy import sqrt, newaxis
from numpy.fft import irfft, rfftfreq
from numpy.random import normal
from numpy import sum as npsum
from scipy import interpolate

# Classe che genera una serie temporale di rumore colorato con la PSD desiderata.
# Parametri :   size =numero di punti della serie temporale che si vuole ottenere
#  ---------  >>>>  log10(asd)= log10(sqrt(PSD)) deve essere un oggetto di tipo interp1d
#               f max e fmin = frequenza massima e minima con cui si interpola la asd
#
#               La classe campiona la asd con la frequenza di campionamento della serie temporale
#                ---->N.B. vedere teo campionamentp
#######################
#-----------------------<--N.B--> si assume che la asd di input sia in realtà log10 della asd "vera"
#-----------------------<--N.B-->  questo per migliorare l'nterpolazione lineare
#-----------------------<--N.B-->  la serie temporale generata ha asd sovrapponibile a quella "vera"
#######################
#!!! Nel caso in cui size=None questa classe assume come size la lunghezza di times array




class GWColoredNoise(GWNoiseSource):

    def __init__(self,name,size,asd,fmin,fmax):
        GWNoiseSource.__init__(self, name, t0=0 ,tau=0)
        self.size = size
        self.asd=asd
        self.fmin=fmin
        self.fmax=fmax

    #(size, asd,fmin, fmax)
    #  !!!!     asd oggetto di tipo "interp1d"
    #  !!!!    max=frequenza massima della psd, per virgo 8180
    def __call__(self, times_array,dect=None):




        num=len(times_array)
        """
        Based on the algorithm in:
        Timmer, J. and Koenig, M.:
        On generating power law noise.
        Astron. Astrophys. 300, 707-710 (1995)
        fmin : float, optional
        Low-frequency cutoff.
        Default: 0 corresponds to original paper. It is not actually
        zero, but 1/samples.
        """
        try:
            size = list(self.size)
        except TypeError:
            size = [self.size]
            # The number of samples in each time series
            samples = size[-1]

        # Use fft functions for real output (-> hermitian spectrum)
        #print("!!!!!!",1 / (2 * self.fmax))
        if self.size is None:
            samples=len(times_array)
            sample_rate=1./(times_array[1]-times_array[0])

        f = rfftfreq((samples), 1. / ( sample_rate) ) # dt calcolato come 1/2f_max e fmax=8100
        #print("ff",f)
        #print("f_max_interp",np.amax(f))
        #print("f_min_interp",np.amin(f))

        # Build scaling factors for all frequencies
        s_scale = np.zeros(samples)
        fmin = max(self.fmin, 1. / samples)  # Low frequency cutoff
        ix = npsum(s_scale < fmin)  # Index of the cutoff
        if ix and ix < len(s_scale):
            s_scale[:ix] = s_scale[ix]
        f = f[f <= self.fmax]
        f = f[f >= self.fmin]
        asd_new = self.asd(f)
        asd_new = 10**asd_new
        df = np.abs(f[1] - f[2])

        s_scale = ((0.5)) * asd_new * df ** 0.5

    # Adjust size to generate one Fourier component per frequency
        size[-1] = len(f)

    # Add empty dimension(s) to broadcast s_scale along last
    # dimension of generated random power + phase (below)
        dims_to_add = len(size) - 1
        s_scale = s_scale[(newaxis,) * dims_to_add + (Ellipsis,)]

    # Generate scaled random power + phase
    # print(s_scale)
        sr = normal(scale=(s_scale), size=size)
        si = normal(scale=(s_scale), size=size)
    # print(sr)
    # If the signal length is even, frequencies +/- 0.5 are equal
    # so the coefficient must be real.
        if not (samples % 2): si[..., -1] = 0

    # Regardless of signal length, the DC component must be real
        si[..., 0] = 0

    # Combine power + corrected phase to Fourier components
        s = sr + 1J * si

    # Transform to real time series
        y = samples * irfft(s, n=samples, axis=-1)
    # print(f)
        #print("len",len(y))
        return y


