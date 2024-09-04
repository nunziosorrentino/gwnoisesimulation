import unittest
import numpy as np
import logging

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation

from gwskysim.utilities.gwtools import apply_delay
from gwskysim.detectors.gwdetector import GWDetector
from gwskysim.sources.gwperiodicsource import GWPeriodicSource

logging.basicConfig(level=logging.INFO)

params = {
            'name': 'basic source',
            'ra': 0,
            'dec': 0,
            'polarization': 2.13,
            'distance': 3.086e+21,
            'inclination': 0.35,
            'ecc': 10e-5,
            'ephemeris': [9.63939485809138, -8.894189277846e-12, 1.0887773549E-21],
            'wob_angle': 0.13,
            'mom_inertia1': 10e45,
            'initial_phase': 1
        }
source = GWPeriodicSource(**params)
detector = GWDetector('V1')

class TestCases(unittest.TestCase):

    def test_gw_travel_time_correction(self):

        # test if fix_elem works properly
        times = np.linspace(0, 1100)
        new_times = apply_delay(times, source, detector, fix_elem=0)
        self.assertEqual(times[0], new_times[0], msg='fix_elem does not work properly (elem 0)')

        times = np.linspace(0, 1400, 360)  # check over a year
        other_new_times = apply_delay(times, source, detector, fix_elem=32)
        self.assertEqual(times[32], other_new_times[32], msg='fix_elem does not work properly (elem 32)')