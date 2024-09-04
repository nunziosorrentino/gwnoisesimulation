import numpy as np
import itertools
from inspect import getmembers, isfunction
#from gwskysim.models.closenc_wavefuncs import shape_func as sf
#from gwskysim.models.closenc_wavefuncs import shape_matrix as sm
from gwskysim.models import closenc_wavefuncs
from gwskysim.sources.gwpointsource import GWPointSource
from gwskysim.utilities.gwsim_util import pc_to_Msun, Msun_to_sec

# Get close encounters wave functions
dict_f = dict(getmembers(closenc_wavefuncs, isfunction))
del dict_f['sintable']
del dict_f['costable']
del dict_f['th_ci_si']
del dict_f['cb_sb']

# Define one function for all
def shape_func(c_or_s, p_or_c, phase, incl, pol):
    # name of the function in global scopes
    def shape_partial(n, k):
        func_str = '{}{}_{}_{}'.format(c_or_s, p_or_c, k, n)
        # if the function doesn't exist, it returns zero
        try:
            return dict_f[func_str](phase, incl, pol)
        except KeyError:
            return np.full(len(phase), 0)
    return shape_partial        

# Define a relative matrix with indexes 'k' and 'n'
def shape_matrix(c_or_s, p_or_c, phase, incl, pol):
    #import time
    #i_t = time.time()
    shape_f = shape_func(c_or_s, p_or_c, phase, incl, pol)
    sm_ = np.array([[shape_f(n, k) for k in range(7)] for n in range(3)])    
    #k_n = list(itertools.product(np.arange(7), np.arange(3)))
    #k_n = np.array(k_n).T
    #sm_ = np.fromiter([[shape_f(n, k) for k in range(7)] for n in range(3)], float, count=21)
    #print('AAAAA', sm_.shape)
    #f_t = time.time()
    #print('CHECK TIME shape matrix: {} s'.format(f_t - i_t))     
    return sm_

class GWBinarySource(GWPointSource):
    """Binary inspiraling eccentric source (BHBH, NSBH and NSNS) 
       that emit as a close encounter.

     Arguments
     ---------
     name : string
        The name of the source.

     ra : float or array-like of floats
        Right ascension of the source in degrees.

     dec : float or array-like of floats
        Declination of the source in degrees. 

     polarization : float or array-like of floats
        Polarization angle of the source emission in degrees. 
   
     inclination : float or array-like of floats
        Inclination of the source main axis in degrees.
        
     distance : float or array-like of floats
        Distance of the source in parserc.   

     e0 : float or array-like of floats
        Eccentricity of the orbit containing the neutron star.

     p_0 : float or array-like of floats
        Semi-latus rectum of the binary system over total mass (G=c=1).

     m1 : float or array-like of floats
        Mass of the first body of the system in solar mass units.

     m2 : float or array-like of floats
        Mass of the second body of the system in solar mass units.

     tp : float or array-like of floats
        Epoch of periastron passage in GPS time in seconds.
    """
    def __init__(self, name, ra, dec, polarization, inclination, distance, e0, p_0, m1, m2, tp=0.):
        """Constructor.
        """
        GWPointSource.__init__(self, name, ra, dec, polarization, inclination)
        distance = pc_to_Msun(distance)
        self.distance = distance
        self.e0 = e0
        self.p_0 = p_0
        self.m1 = m1
        self.m2 = m2 
        # Time of periatron passage
        self.tp = tp   

        # total mass, reduced mass and ... mass
        M = m1 + m2 
        mu = m1*m2/M
        eta = mu/M

        # reduced semi-latus rectum
        p0 = p_0*M
        self.p0 = p0

        # high eccentricity parameter
        eps0 = 1 - (e0**2)

        # usefull quantities
        nblock_1 = 1/M
        nblock_2 = (eps0/p_0)**(3/2)
        n0 = nblock_1*nblock_2
        
        fblock_1 = 96/(10*np.pi) 
        fblock_2 = eta/(M*(p_0**4))
        fblock_3 = np.sqrt(eps0)
        fblock_4 = 1 + (73/24)*(e0**2) + (37/96)*(e0**4)
        Frr = fblock_1*fblock_2*fblock_3*fblock_4
    
        # define related attibutes 
        self.M = M
        self.mu = mu
        self.eta = eta
        self.eps0 = eps0
        self.n0 = 1/Msun_to_sec(1/n0)
        #print('CHECK this, Torb = ', Msun_to_sec((2*np.pi)/self.n0))
        self.Frr = 1/Msun_to_sec(1/Frr)   

    def t_edge(self, l_=np.pi):
        r"""Considering a single periastron passage, returns epoch referred to
           a certain eccentric anomaly value.

        Note
        ----
        The default eccentric anomaly value is $\pi$, 
        that corrispionds to half a lap.
        """
        assert l_!=0
        #print('aaa', self.Frr, self.n0)
        
        lblock_1 = np.log(((2*np.pi*l_*self.Frr)/self.n0) + 1) 
        if np.isnan(lblock_1):
            if l_<0:
                return self.tp - 30 
            if l_>0:
                return self.tp + 30             
        lblock_2 = 2*np.pi*self.Frr
        #print('bbb', (lblock_1/lblock_2) + self.tp)
        return (lblock_1/lblock_2) + self.tp

    # Set initial time in format you prefer 
    # l is eccentric anomaly
    def ecc_anomaly(self, times_array):
        """
        Eccentric anomaly of the first body at each time.
        """
        dt = np.subtract(times_array, self.tp)
        block1 = self.n0/(2*np.pi*self.Frr)      
        block2 = np.subtract(np.exp(2*np.pi*self.Frr*dt), 1)     
        return block1*block2

    # How orbital parameters change with time during pericenter passage
    def ecc(self, times_array):
        """
        Eccentricity of the orbit at each time.
        """
        block1 = 20.267*self.eta*self.e0/np.power(self.p_0, 2.5)
        block2 = 1 + 0.39803*np.power(self.e0, 2)
        return self.e0 - block1*block2*self.ecc_anomaly(times_array)

    def eps(self, times_array):
        """
        High eccentricity parameter at each time.
        
        eps = 1 - ecc**2
        """
        return np.subtract(1, np.power(self.ecc(times_array), 2))

    def semi_lat_rect(self, times_array):
        """
        Semi latus rectum of the orbit at each time.
        """
        block1 = 12.8*self.eta/np.power(self.p_0, 2.5)
        block2 = 1 + 0.875*np.power(self.e0, 2)
        return self.M*self.p_0*(1 - block1*block2*self.ecc_anomaly(times_array))

    def phase(self, times_array):
        """
        Phase of the first body at each time.
        """
        sqrt_epst = np.sqrt(self.eps(times_array))
        return self.ecc_anomaly(times_array)/(np.log((1+sqrt_epst)/self.ecc(times_array))-sqrt_epst)

    def ch(self, k, times_array):
        """
        Composed hyperbolic function that defines the scale of binary waveform functions (see N. Loutrel 2019).
        """
        return np.cosh((k/3)*np.arcsinh(self.phase(times_array)))

    def sh(self, k, times_array):
        """
        Composed hyperbolic function that defines the scale of binary waveform functions (see N. Loutrel 2019). 
        """
        return np.sinh((k/3)*np.arcsinh(self.phase(times_array)))
    
    #@profile
    def __call__(self, times_array):
        """
        This call method takes the GW template of a close-encounters system defined in models module
        given a time array as input.
        """
        #import time
        #i_time = time.time()   
        tcut_max = self.tp + 0.5
        tcut_min = self.tp - 0.5
        if self.t_edge(np.pi/4) - self.t_edge(-np.pi/4) < 1.:
            tcut_max = self.t_edge(np.pi/4)
            tcut_min = self.t_edge(-np.pi/4)        
        times_array = times_array[times_array <= tcut_max]
        times_array = times_array[times_array >= tcut_min]
        if np.array(times_array).size==0:
            return None    	
        # the units are with G=c=1
        # Now the amplitude function
        ampl = -(self.eta*(self.M**2))/(self.semi_lat_rect(times_array)*self.distance)
        # Now the same for the scale functions
        scale_ = np.array([self.eps(times_array)**n for n in range(3)])
        # Now the hyperbolic functions
        hype_c = self.ch(np.arange(7), np.array(np.split(np.tile(times_array, 7), 7)).T)
        hype_c = hype_c.T 
        hype_s = self.sh(np.arange(7), np.array(np.split(np.tile(times_array, 7), 7)).T)
        hype_s = hype_s.T
        #f_time = time.time()
        #print('CHECK TIME first part of call: {} s'.format(f_time-i_time))  
        #i_time = time.time() 
        # Define plus polarization
        inclin_ = np.deg2rad(self.inclination)
        polar_ = np.deg2rad(self.polarization)
        h_p = shape_matrix('c', 'p', self.phase(times_array), inclin_, polar_)*hype_c
        h_p += shape_matrix('s', 'p', self.phase(times_array), inclin_, polar_)*hype_s
        h_p = np.swapaxes(np.swapaxes(h_p, 0, 1)*scale_, 0, 1)              
        h_p = ampl*h_p.sum(axis=(0, 1))  
        # Define cross polarization              
        h_c = shape_matrix('c', 'c', self.phase(times_array), inclin_, polar_)*hype_c
        h_c += shape_matrix('s', 'c', self.phase(times_array), inclin_, polar_)*hype_s
        h_c = np.swapaxes(np.swapaxes(h_c, 0, 1)*scale_, 0, 1)               
        h_c = ampl*h_c.sum(axis=(0, 1))        
        #f_time = time.time()
        #print('CHECK TIME h_p, h_c: {} s'.format(f_time-i_time)) 
        #print('bbb', h_p, h_c, times_array) 
        return h_p, h_c, times_array[0]      
        
if __name__ == '__main__':

    import time
    import matplotlib.pyplot as plt
    from scipy.fft import fft, fftfreq
    from gwskysim.sources.gwwhitenoise import White_gaussian_noise
    m1_v = 1.4
    m2_v = 1.4
    dist_v = 5e8 #in pc
    p0_v = 25 # in M unit
    e0, p_0, m1, m2, dist, incl, pol = (0.95, p0_v, m1_v, m2_v, dist_v, 0, 0)
    t0_p = 1241921261.593
    gw_source = GWBinarySource('CE {}-{} Msun'.format(m1, m2), 295.851965, 
                               -15.5165, pol, incl, dist, e0, p_0, m1, m2, tp=t0_p)
    tl_ = gw_source.t_edge(np.pi/4)
    tl_m = gw_source.t_edge(-np.pi/4)
    print("Signal duration: {} s".format(tl_-tl_m))
    t_f = tl_m 
    duration = tl_- tl_m
    #print('aaaa', duration, tl_, tl_m)
    fs = 4096
    
    init_time_1 = time.time()
    hp, hc = gw_source._import_template(np.linspace(t_f, t_f+duration, int(duration*fs)))
    end_time_1 = time.time()
    print('"_import_template" RUN TIME {} s.'.format(end_time_1 - init_time_1))

    init_time = time.time()
    h_1 = gw_source.generate_signal(t_f, duration, fs, detectors = ['V1', 'L1', 'H1'])
    end_time = time.time()
    print('"gen_signal" RUN TIME {} s.'.format(end_time - init_time))
    
    plt.figure()   
    plt.plot(hp.get_sampled_time()-t0_p, hp.sampled_strain/hp.sampled_strain.max(), label=r'$h_{+}$')
    plt.title(r'CE {}-{} Msun p={} Msun e={}'.format(m1, m2, p_0, e0))
    plt.ylabel('Nornalized Strain []')
    plt.xlabel('Time from {} GPS [s]'.format(t0_p)) 
    plt.legend()
    plt.tight_layout()   
    
    ffts_ = fft(hp.sampled_strain)
    N = int(duration*fs)
    print('aaa', N//2)
    T = 1.0/fs
    f_reqs = fftfreq(N, T)[:N//2]
    print('bbb', len(f_reqs))
    amplitudes = 2.0/N * np.abs(ffts_[0:N//2])
    print('ccc', len(amplitudes))
    plt.figure()
    plt.title(r'FFT CE {}-{} Msun p={} Msun e={}'.format(m1, m2, p_0, e0))
    plt.plot(f_reqs[:len(amplitudes)], amplitudes/amplitudes.max(), label=r'$h_{+}$') 
    plt.xlim((30, 1200)) 
    plt.ylim((0, 0.6))
    plt.legend()
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Normalized strain []')
    #h_1[0].plot()
    #h_1[1].plot()
    #h_1[2].plot()
    plt.show()
        
