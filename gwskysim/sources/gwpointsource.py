import numpy as np
from gwskysim.sources.gwsource import GWSource
from gwskysim.detectors.gwdetector import GWDetector

from gwskysim.sources.gwsignal import GWSignal
from gwskysim.detectors.gwdetector import GWDetector
from gwskysim.utilities.gwtools import apply_delay

class GWPointSource(GWSource):
    """General point Gw source located in the sky.

    Note
    ----
    ra, dec and polarization are sky position of GW 
    and polarization is the polarization angle of tensor wave (Thorn 1987)

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
   
     inclination : float or array-like of floats
        Inclination of the source main axis in radians. 
    """
    def __init__(self, name, ra, dec, polarization, inclination, **kwargs):
        """Constructor.
        """
        GWSource.__init__(self, name)
        self.ra = ra
        self.dec = dec
        self.polarization = polarization
        self.inclination = inclination
        self.__dict__.update(kwargs)   

    def generate_signal(self, m_start_time, m_duration, m_sample_rate,
                        detectors=None, approximant=None, lower_freq=40):
        """
        Generate the strain signal in the proper detector using the
        specified __call__ method of the derivate class.

        :param m_start_time:
        :param m_duration:
        :param m_sample_rate:
        :param detectors: list of strings
            Only string lists are allowed, i.e.: ['V1', 'H1']
        :param approximant:
        :param lower_freq:
        :return:
        """

        if detectors is None:
            detectors = ['V1', 'L1', 'H1']

        # set times array
        gpsts = np.linspace(m_start_time,
                        m_start_time+m_duration,
                        int(m_duration*m_sample_rate))

        # initialize a list of h for every given detector
        hl = []
        for detname in detectors:

            # initialize detector instance
            det = GWDetector(detname)

            # times corrected with delay
            #ct = gpsts
            ct = apply_delay(gpsts, self, det, input_scale='utc', fix_elem=0)                              

            # import template section
            hp, hc = self._import_template(ct, approximant,
                                           lower_freq, m_sample_rate)
            if (hp, hp)==(None, None):
                return None                               

            # projection
            signal_strain = det.project_wave(hp, hc, self.ra, self.dec,
                                              self.polarization)

            hl.append(signal_strain)

        return hl

