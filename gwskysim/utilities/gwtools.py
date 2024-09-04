import numpy as np
import logging

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from scipy.interpolate import interp1d


def apply_delay(times, source, detector, input_scale='tai', fix_elem=None, max_delta_t=1200):
    """
    Adjust the time of arrival of the source signal taking into account the travel time from the detector to the
    SSB, including relativistic corrections. Delay is evaluated once every (at most) max_delta_t seconds;
    the other points are inferred using interpolation.

    :param times: array of float
        The time(s) as measured in the detector site (in 'gps' format).
    :param source: GWPointSource object
        The source object
    :param detector: GWDetector object
        The detector object
    :param input_scale: string
        The input scale can be chosen accordingly with the list of available scales in Astropy documentation. By default,
        the scale is 'tai'.
    :param fix_elem: int or None
        Time array index which points the element that the user wants to remain unchanged. This feature applies (in
        general) a uniform shift on the whole corrected time array. Default value is None (no shift applied).
    :param max_delta_t: float
        Delay will evaluated at most every max_delta_t seconds. Default value is 1200 (sec): this parameter can be
        reduced, but not incremented as it can give rise to incorrect results.
    :return: array of float
        The times array in SSB in 'tdb' scale.
    """

    assert max_delta_t <= 1200, 'It is not allowed to set max_delta_t above 1200 sec. The distance between delay ' \
                               'samples is too high.'

    # create astropy SkyCoord object and SkyCoord object
    source_coords = SkyCoord(source.ra, source.dec, frame='icrs', unit='deg')
    det_location = EarthLocation.from_geodetic(detector.longitude, detector.latitude)

    # the piece of code below ensure that interpolation has at least 10 points
    tmp = 10
    if times[-1] - times[0] < tmp * max_delta_t:
        num_of_delay_samples = tmp + 1  # plus one for end point
    else:
        num_of_delay_samples = int((times[-1] - times[0]) / max_delta_t) + 2  # start and end points
    times_s = np.linspace(times[0], times[-1], num_of_delay_samples)
    assert times_s[1] - times_s[0] < max_delta_t

    apply_delay.num_of_delay_samples = num_of_delay_samples  # set attribute

    # create astropy Time object for sample points and the relative delay, then correct times
    a_times_s = Time(times_s, format='gps', scale=input_scale)
    a_delays_s = a_times_s.light_travel_time(source_coords, location=det_location)
    correct_times_s = (a_times_s.tdb + a_delays_s).value

    # interpolate to find correct times
    delay_func = interp1d(times_s, correct_times_s)
    correct_times = delay_func(times)

    if fix_elem is not None:
        shift = correct_times[fix_elem] - times[fix_elem]
        correct_times -= shift

    return correct_times


if __name__ == '__main__':

    import logging
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.time import Time
    import time

    from gwskysim.detectors.gwdetector import GWDetector
    from gwskysim.sources.gwperiodicsource import GWPeriodicSource

    logging.basicConfig(level=logging.INFO)

    # create source and detector
    params = {
        'name': 'fake_source',
        'ra': 0,
        'dec': 0,
        'polarization': 2.13,
        'distance': 3.086e+21,
        'inclination': 0.35,
        'ecc': 10e-5,
        'ephemeris': [1],
        'wob_angle': 0.13,
        'mom_inertia1': 10e45,
        'initial_phase': 1
    }
    source = GWPeriodicSource(**params)
    detector = GWDetector('V1')

    # initialize times
    # times = np.linspace(0, 31536000, 1000000)  # check over a year
    times = np.linspace(0, 172800, 3600)  # check over a day
    times = times + 1268740818  # start with march equinox 2020

    # evaluate delay
    t0 = time.time()
    new_times = apply_delay(times, source, detector)
    t1 = time.time()
    print('run time: {}'.format(t1-t0))
    delay = new_times - times

    # plot
    plt.figure()
    plt.title('Delay')
    plt.xlabel('seconds from starting time $(s)$')
    plt.ylabel('delay (s)')
    plt.plot(times - times[0], delay)
    plt.show()

    # try different duration
    a = 100
    b = 10000 + np.arange(100) * a

    num_of_delay_samples_list = []
    for duration in b:
        times = np.linspace(0, duration)
        new_times = apply_delay(times, source, detector)
        num_of_delay_samples_list.append(apply_delay.num_of_delay_samples)

    plt.figure()
    plt.plot(b, num_of_delay_samples_list, '.')
    plt.title('delay sample points')
    plt.xlabel('duration $(s)$')
    plt.ylabel('points')
    plt.show()