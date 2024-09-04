
Downloading data from GWOSC
==========================

In order to download real data from the GWOSC, enter the ``/some/path/gwskysim/config``
directory and open the "download_raw_data_config.json" file. This is a JSON configuration
file which contains the following keys:

``observation_run``
    The observation run of which you want to download the data (and the sampling rate).
    Default is "O2_4KHZ_R1" which corresponds to O2 sampled at 4192 Hz. Check
    ``www.gw-openscience.org`` for other possibilities.

``detectors``
    The detectors whose data you want to download. Default is "[H1,L1,V1]".

``GPS_start_time`` and ``GPS_end_time``
    The GPS time interval in which you want to download the data.

``min_duty_cycle``
    Only files with (at least) this duty cycle will be downloaded. Must be a number
    in range [0,100]. Default is 50.

.. note::
    Only data where ALL the specified detectors were working with the minimum duty cycle
    will be downloaded.

Once you've adjusted the JSON configuration file, return to your work directory
``/some/path/`` (the directory in which gwskysim is contained) and run the following
command:

.. code-block:: bash

   download_raw_data.py ./gwskysim/config/download_raw_data_config.json

Data will be downloaded and stored in a directory called with the standardized name
"raw_data_observationrun_starttime_endtime". Data will be divided in sub-directories,
one for each detector.

.. note::
    GWOSC data can take up a lot of memory (several GB). Use debug mode (-d) if you
    want to check the list of files that will be downloaded before actually download them.

Select noise according to quality
==========================

Once you have downloaded the data you can do a quality check to create a list of "good"
noise GPS times that will allow you to quickly download the detectors noise (e.g. through
a python generator). In order to do that, enter the ``/some/path/gwskysim/config``
directory and open the "generate_noise_times_config.json" file.
This is a JSON configuration file which contains the following keys:

``noise_directory``
    The directory that contains the data downloaded from the GWOSC (see above), should be
    of the form "raw_data_observationrun_starttime_endtime".

``detectors``
    The detectors for which you want to check the quality. Default is "[H1,L1,V1]".

``output_directory_name``
    The directory in which you want to create the `*.hdf5` output file. You can use the
    "noise_directory".

``len_noise_seconds``
    Must be an integer. Noise segments with the chosen quality requirements must have
    at least this length in seconds. Default is 128. Segments of shorter length are
    discarded.

``data_quality_bits``
    The Data Quality Bits. You can find specific information on GWOSC website.
    The more bits you select, the better your data quality will be. Default is
    "data_quality_bits": [0, 1, 2, 3, 4, 5, 6] (best possible quality)

``injection_bits``
    The Injection Bits. You can find specific information on GWOSC website.
    Default is "injection_bits": [0, 1, 2, 4]. This means that the only type
    of hardware injection allowed is continuous wave injections.

Once you've adjusted the JSON configuration file, return to your work directory
``/some/path/`` (the directory in which gwskysim is contained) and run the following
command:

.. code-block:: bash

   generate_noise_times.py ./gwskysim/config/generate_noise_times_config.json

This will create an `*.hdf5` file with a list of integer GPS times.
Suppose that a GPS time T is in the list. This means that in the time interval
[T,T+"len_noise_seconds"] the data have the specified quality requirements
in ALL detectors indicated in the configuration file.
