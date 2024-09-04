import numpy as np
import logging
from scipy.special import factorial

from gwskysim.sources.gwpointsource import GWPointSource


class GWPeriodicSource(GWPointSource):
    """A non axisimmetric rotating neutron star.

     Arguments
     ---------
     name : string
        The name of the source.

     ra : float or array-like of floats
        Right ascension of the source in degrees.

     dec : float or array-like of floats
        Declination of the source in degrees. 

     polarization : float or array-like of floats
        Polarization angle of the source emission in radians. 

     distance : float or array-like of floats
        Distance of the source in solar mass units (G=c=1).
   
     inclination : float or array-like of floats
        Inclination of the source main axis in radians.

     ecc : float or array-like of floats
        Eccentricity of the neutron star body.
        
     ephemeris : array-like of floats
        The parameters are interpreted as follows: if a single value is given ([f0]), it is the rotational
        frequency of the source. If multiple values are given [f_0, f_1, ..] then f_0 is the rotational frequency of
        the source, f_1 is the first spindown parameter (fdot), and so on.

     wob_angle : float or array-like of floats
        Angle representing the asimmetry of the neutron star in radians.
        
     mom_inertia1 : float or array-like of floats
        First component of the moment of inertia of the neutron star.

     initial_phase : float or array-like of floats
        Initial phase of the signal in radians.
    """

    def __init__(self, name, ra, dec, polarization, distance, inclination, ecc, ephemeris,
                 wob_angle, mom_inertia1, initial_phase):
        """Constructor.
        """
        GWPointSource.__init__(self, name, ra, dec, polarization, inclination)
        self.ecc = ecc
        self.wob_angle = wob_angle
        self.mom_inertia1 = mom_inertia1
        self.distance = distance
        self.initial_phase = initial_phase
        self.inclination = inclination
        self.ephemeris = np.array(ephemeris)


    def frequency(self, times_array):
        """
        This method evaluates the spin frequency of the source, for every instant of time.

        :param times_array: array which contains the time samples in the selected time range
        :return: frequency array, corresponding to the time sequence
        """

        s = len(self.ephemeris) # number of spindown parameters
        t_length = len(times_array) # number of time samples

        f = np.zeros(t_length) # initialize f array
        for i in range(s):
            f += self.ephemeris[i] * (times_array ** i) / factorial(i)

        return f


    def amplitude(self, times_array):

        """
        A method to evaluate the amplitude of the signal, for every instant of time, in the SSB reference system. It is
        assumed that all the parameters are expressed in CGS units.

        :param times_array: array which contains the time samples in the selected time range
        :return: amplitude array
        """

        # see Jaranowski
        d_0 = (self.ecc / 10e-5) * (self.mom_inertia1 / 10e45) * (3.086e+21 / self.distance)

        return 4.23 * 10e-25 * d_0 * (self.frequency(times_array) / 100) ** 2


    def phase(self, times_array):

        """
        A method to evaluate the phase of the gravitational wave in the SSB reference system, for every instant of time,
        taking into account the spindown parameters.

        :param times_array: numpy ndarray of floats
            array which contains the time samples in the selected time range
        :return: phase array
        """

        # is sampling theorem condition satisfied?
        max_delta_t = np.max(np.diff(times_array))
        tmp = np.max(self.frequency(times_array))
        max_phase_derivative = 4 * np.pi * tmp
        max_delta_phase = max_delta_t * max_phase_derivative

        if max_delta_phase > np.pi: # which means less than two samples for each period
            logging.warning('Sample rate does not exceed the Nyquist rate for the given source: '
                            'this may lead to problems in case one wants to reconstruct the phase in the frequency domain.\n'
                            'Recommended sampling frequency: at least {} Hz'.format(4 * tmp))

        s = len(self.ephemeris)  # number of spindown parameters
        t_length = len(times_array)  # number of time samples

        p = np.full(t_length, self.initial_phase, dtype='float64') # initialize phase array with initial phase
        for i in range(s):
            p += 2 * np.pi * self.ephemeris[i] * (times_array ** (i+1)) / factorial(i+1)

        return 2 * p # GW phase is twice the rotational phase

    def __call__(self, times_array):
        t = times_array
        i = self.inclination

        # see Woan (2005)
        h_p = (1/2) * self.amplitude(t) * (1 + np.cos(i) ** 2) * np.cos(self.phase(t))
        h_c = self.amplitude(t) * np.cos(i) * np.sin(self.phase(t))

        return h_p, h_c


if __name__ == '__main__':
    import time

    print('Creating source...')
    # parameters: random (but typical) values
    params = {
        'name': 'PSR J2021+3651',
        'ra': 305.27275,
        'dec': 36.851333,
        'polarization': 2.13,
        'distance': 3.086e+21,
        'inclination': 0.35,
        'ecc': 10e-5,
        'ephemeris': [9.63939485809138, -8.894189277846e-12, 1.0887773549E-21],
        # 'ephemeris': [2],
        'wob_angle': 0.13,
        'mom_inertia1': 10e45,
        'initial_phase': 1
    }
    source = GWPeriodicSource(**params)

    # initialize detector
    det_name = 'V1'

    # # try
    # times_array = np.linspace(0,10)
    # print('times_array:\n{}'.format(times_array))
    #
    # # test frequency
    # _t0 = time.time()
    # print('Frequency: \n{}'.format(source.frequency(times_array)))
    # _t1 = time.time()
    # print('run time: {}'.format(_t1 - _t0))
    #
    # # test amplitude
    # _t0 = time.time()
    # print('Amplitude: \n{}'.format(source.amplitude(times_array)))
    # _t1 = time.time()
    # print('run time: {}'.format(_t1 - _t0))
    #
    # # test phase
    # _t0 = time.time()
    # print('Phase: \n{}'.format(source.phase(times_array)))
    # _t1 = time.time()
    # print('run time: {}'.format(_t1 - _t0))

    # time range parameters
    t_f = 0
    duration = 360
    fs = 100

    print('Generating signal...')
    _t0 = time.time()
    h_v1 = source.generate_signal(t_f, duration, fs, detectors=['V1'])[0]
    _t1 = time.time()
    print('run time: {}'.format(_t1 - _t0))

    import matplotlib.pyplot as plt

    plt.figure(figsize=(20, 4))
    plt.plot(h_v1.get_sampled_time(), h_v1)
    plt.show()
