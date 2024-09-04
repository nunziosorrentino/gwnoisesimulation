Source Simulation (BBH)
***********************

Now, from the sources parameters we can generate the waveforms.

Enter the ``/some/path/gwskysim/config`` directory and open the
"generate_waveforms_config.json" file. This JSON configuration file includes the
following keys:

``data_directory``
    The directory in which sources parameters generated before are contained. Default is
    the directory "sources_parameters".

``output_directory_name``
    The name of the directory in which waveforms output files will be created. This
    directory will be inside the population directory at the same level as
    "sources_parameters" directory.

``parallel_processes``
    Number of parallel processed that will be used in the simulation. This speeds up
    the simulation.
    The number of processes should be less than the number of cores on your computer.

``approximant``
    The model that will be used to generate the waveforms. Default is "SEOBNRv4".

``sample_freq``
    The sampling frequency (measured in Hz) of the waveforms. Default is "4096".

``lower_freq``
    The starting frequency (in Hz) of each waveform. Default is "16".

``detectors``
    The network of the detector in which you want to simulate the waveforms.
    Default is "[H1,L1,V1]" (Hanford, Livingston, Virgo).

``reference_GPS_time``
    The GPS time at which the signals arrive at geocenter. Default is "1238112018.0"
    that is the starting time of O3.

Once you've adjusted the JSON configuration file, go inside the population directory
``/some/path/population_directory`` (same level as "sources_parameters") and run
the following command:

.. code-block:: bash

   generate_waveforms.py ../gwskysim/config/generate_waveforms_config.json

At the end of the script you will see the message:

.. code-block:: console

    All waveforms have been generated and saved to disk

.. note::
     This could take many hours to complete!

The script will create a directory called "waveforms". The output `*.hdf5` files
will be contained in sub-directories, one for each detector used in the simulation.
By default you will find three directories: "H1", "L1", "V1".
For a given source, you will have three waveforms, one for each detector, projected
in the detector antenna pattern. The time delay due to distance between
detectors is taken into account.

