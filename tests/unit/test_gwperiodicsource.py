import unittest
import numpy as np

from gwskysim.sources.gwperiodicsource import GWPeriodicSource


class TestCases(unittest.TestCase):

    def test_amplitude(self):

        params = {
            'name': 'PSR J2021+3651',
            'ra': 305.27275,
            'dec': 36.851333,
            'polarization': 2.13,
            'distance': 3.086e+21,
            'inclination': 0.35,
            'ecc': 10e-5,
            'ephemeris': [1],
            'wob_angle': 0.13,
            'mom_inertia1': 10e45,
            'initial_phase': 1
        }
        basic_source = GWPeriodicSource(**params)
        times_array = np.linspace(0, 3600, 36000)

        # amplitude does not vary if f is constant
        amplitude = basic_source.amplitude(times_array)
        self.assertEqual(amplitude[0], amplitude[32], msg='amplitude varies unexpectedly')
