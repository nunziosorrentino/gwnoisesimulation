import numpy as np
from astropy.time import Time
from gwskysim.sources.gwsignal import GWSignal

class GWDetector:
    """
    Class for a gravitational wave detector

    Arguments
    ---------
    name : string
        Name of the detector. Use GWDetector.get_detector_name to see a list of available values.
    """

    # in order, we have:
    # [latitude, longitude, orientation of detector arms, angle between the arms]
    available_detectors = {
        'V1': [43.63, -10.5, 116.5, 90.],
        'H1': [46.45, 119.41, 171.8, 90.],
        'L1': [30.56, 90.77, 243., 90.],
        'GEO600': [52.25, -9.81, 68.775, 94.33],
        'TAMA300': [35.68, -139.54, 225., 90.],
        'ET': [40.44, -9.4566, 116.5, 60.] # Sardinia site hypothesis
    }

    def __init__(self, name, t_ref=0):
        """Constructor"""

        self.name = name
        self.t_ref = t_ref

        if name not in self.available_detectors.keys():
            raise ValueError("Not valid argument ({}) for 'name' parameter.".format(name))

        self.latitude = self.available_detectors[name][0]
        self.longitude = self.available_detectors[name][1]
        self.gamma = self.available_detectors[name][2]
        self.zeta = self.available_detectors[name][3]

    @staticmethod
    def _ab_factors(g_, lat, ra, dec, lst):
        """
        Method that calculates the amplitude factors of plus and cross
        polarization in the wave projection on the detector. 
        :param g_: float
            this represent the orientation of the detector's arms with respect to local geographical direction, in
            rad. It is measured counterclock-wise from East to the bisector of the interferometer arms.
        :param lat: float
            longitude of the detector in rad.
        :param ra: float
            Right ascension of the source in rad.
        :param dec: float
            Declination of the source in rad.
        :param lst: float or ndarray
            Local sidereal time(s) in rad. 
        :return: tuple of float or np.ndarray
            relative amplitudes of hplus and hcross.       
        """
        a_ = (1/16)*np.sin(2*g_)*(3-np.cos(2*lat))*(3-np.cos(2*dec))*np.cos(2*(ra - lst))-\
             (1/4)*np.cos(2*g_)*np.sin(lat)*(3-np.cos(2*dec))*np.sin(2*(ra - lst))+\
             (1/4)*np.sin(2*g_)*np.sin(2*lat)*np.sin(2*dec)*np.cos(ra - lst)-\
             (1/2)*np.cos(2*g_)*np.cos(lat)*np.sin(2*dec)*np.sin(ra - lst)+\
             (3/4)*np.sin(2*g_)*(np.cos(lat)**2)*(np.cos(dec)**2)

        b_ = np.cos(2*g_)*np.sin(lat)*np.sin(dec)*np.cos(2*(ra - lst))+\
             (1/4)*np.sin(2*g_)*(3-np.cos(2*lat))*np.sin(dec)*np.sin(2*(ra - lst))+\
             np.cos(2*g_)*np.cos(lat)*np.cos(dec)*np.cos(ra - lst)+\
             (1/2)*np.sin(2*g_)*np.sin(2*lat)*np.cos(dec)*np.sin(ra - lst)

        return a_, b_

    @classmethod
    def get_detector_names(cls):
        return cls.available_detectors.keys()

    @classmethod
    def new_detector(cls, name, latitude, longitude, gamma, zeta):
        """
        This method provides a tool to create a new detector given the parameters.

        :param name: string
            the name of the detector.
        :param latitude: float
            latitude of the detector in degrees.
        :param longitude: float
            longitude of the detector in degrees.
        :param gamma: float
            this represent the orientation of the detector's arms with respect to local geographical direction, in
            degrees. It is measured counterclock-wise from East to the bisector of the interferometer arms.
        :param zeta: float
            angle between the interferometer arms, in degrees.
        :return: GWDetector object
            the detector object.
        """
        assert 90 >= latitude >= -90, 'Not valid value for latitude. Accepted values: -90 < latitude < 90'
        assert 180 >= longitude >= -180, 'Not valid value for longitude. Accepted values: -180 < longitude < 180'
        assert 0 <= gamma <= 360, 'Not valid value for orientation of the detector arms. ' \
                                'Accepted values: 0 < gamma < 360'
        assert 0 <= zeta <= 360, 'Not valid value for the angle between the detector arms. ' \
                               'Accepted values: 0 < zeta < 360'
        cls.available_detectors[name] = [latitude, longitude, gamma, zeta]
        return cls(name)

    def lst_estimate(self, GPS_time):
        """
        Estimate the local apparent sidereal time (in deg) 
        at the observatory site. This is necessary fit the detector arms
        relative rotation over the observation time.
        
        :param GPS_time: float, int, list or np.ndarray 
            time of arrival of the source signal.
        :return: float or np.ndarray
            local sidereal time(s) in rad
        """
        t_ = Time(np.array(GPS_time)-self.t_ref, format='gps', scale='utc',
                  location=('{}d'.format(self.longitude),
                            '{}d'.format(self.latitude)))
        lst_ = t_.sidereal_time('mean').rad    
        return lst_

    def antenna_pattern_functions(self, right_ascension, declination, polarization, GPS_time):
        '''
        Evaluate the antenna pattern functions.

        :param right_ascension: float
            Right ascension of the source in degree.

        :param declination: float
            Declination of the source in degree.

        :param polarization: float
            Polarization angle of the wave in degree.

        :param GPS_time: float, int, list or np.ndarray 
            time of arrival of the source signal.
            
        :return: tuple of float or np.ndarray
            fplus and fcross. 
        '''
        ra = np.deg2rad(right_ascension)
        dec = np.deg2rad(declination)
        pol = np.deg2rad(polarization)
        lat = np.deg2rad(self.latitude)
        g_ = np.deg2rad(self.gamma)
        z_ = np.deg2rad(self.zeta)
        #lst = np.deg2rad(self.lst_estimate(GPS_time))
        lst = self.lst_estimate(GPS_time)

        ampl11, ampl12 = self._ab_factors(g_, lat, ra, dec, lst)

        if self.name == 'ET':
            ampl21, ampl22 = self._ab_factors(g_+(2/3)*np.pi, lat, ra, dec, lst)
            ampl31, ampl32 = self._ab_factors(g_+(4/3)*np.pi, lat, ra, dec, lst)

        fplus = np.sin(z_)*(ampl11*np.cos(2*pol) + ampl12*np.sin(2*pol))
        fcross = np.sin(z_)*(ampl12*np.cos(2*pol) - ampl11*np.sin(2*pol))

        if self.name == 'ET':
            fplus2 = np.sin(z_)*(ampl21*np.cos(2*pol) + ampl22*np.sin(2*pol))
            fplus3 = np.sin(z_)*(ampl31*np.cos(2*pol) + ampl32*np.sin(2*pol))
            fplus = np.sqrt(fplus**2 + fplus2**2 + fplus3**2)

            fcross2 = np.sin(z_)*(ampl22*np.cos(2*pol) - ampl21*np.sin(2*pol))
            fcross3 = np.sin(z_)*(ampl32*np.cos(2*pol) - ampl31*np.sin(2*pol))
            fcross = np.sqrt(fcross**2 + fcross2**2 + fcross3**2)

        return fplus, fcross

    def project_wave(self, hp, hc, ra, dec, pol):
        """
        Calculate the signal stain projected on the detector.
        
        :param hp: GWSignal object
            Plus polarization of the wave in GWSignal format.

        :param hc: GWSignal object
            Cross polarization of the wave in GWSignal format.

        :param ra: float
            Right ascension of the source in deg.

        :param dec: float
            Declination of the source in deg.

        :param pol: float
            Polarization angle of the wave in deg.
            
        :return: GWSignal object
            Final gravitational wave signal (detector stain).          
        """
        assert(len(hp)==len(hc))
        GPS_times = hp.get_sampled_time()
        f_plus, f_cross = self.antenna_pattern_functions(ra, dec, pol, GPS_times)
        
        h_signal = hp*f_plus + hc*f_cross
        h_signal.detector = self.name
        
        return h_signal
        
