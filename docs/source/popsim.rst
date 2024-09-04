
Simulating a population of BBH
==========================


Once you have set up the Python environment you can start generating your own
synthetic gravitational-wave signals!
(at the moment only populations of binary black holes (BBH) are available)

To create a population of BBH, enter the ``/some/path/gwskysim/config`` directory and
open the "bbh-pop-config1.json" file. This is a JSON configuration file that includes
options for the range of all population parameters that you can adjust according
to your needs:

``output_dir`` and ``population_name``
    "population_name" is the name you wish to give to your population.
    "output_dir" is the name of the **population directory**. If you write "null"
    (without the quotation marks) a directory with the standard name
    "GWsim_population_name_yyyymmdd" will be created.
    The output files will be created in a sub-directory of the population directory
    called "sources parameters".

.. note::
     If you have already created a population and you want to increase the number
     of sources of the population (without creating a new one) you can simply write
     the name of the existing population directory as "output dir" parameter and
     follow the normal steps (specified below). This will also generate waveforms for
     the new sources only when running the source simulation.

``random_seed``
    The seed for the random number generator which is needed for a matter of
    reproducibility.

  .. note::
     Be careful, and change the seed when generating ***different*** populations!

``source_type``
    At the moment only "BBH" is available, sorry =(

``sources_per_file`` and ``n_files``
    The output will be organized in a certain number of files (parameter "n_files") and
    each file will have a "sources_per_file" number of sources. This organization is
    essential in order to avoid memory issues when using the data for machine learning
    purposes.

``parameters``
    Distribution and Range (min and max value) for each parameter. Default distribution
    is "uniform" for all parameters except "inclination" which is cosine-distributed
    ("cos") and declination which is sin-distributed ("sin").

Once you've adjusted the JSON configuration file, return to your work directory
``/some/path/`` (the directory in which gwskysim is contained) and run the following
command:

.. code-block:: bash

   generate_sources_parameters.py ./gwskysim/config/bbh-pop-config1.json

At the end of the script you will see the message:

.. code-block:: console

  Sources parameters files have been successfully created

This script creates the population directory (what we called "output_dir"). The output
`*.hdf5` files will be contained in a sub-directory of the population directory called
"sources_parameters".
Each file ("n_files" in total) contains a table with a certain number of sources
("n_sources_per_file") randomly sampled from the parameter space.

