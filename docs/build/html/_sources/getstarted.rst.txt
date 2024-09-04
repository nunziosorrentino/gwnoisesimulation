Getting Started
*******************
Before starting the simulations, you have to initialize and setup the package

Initialization
===================

Before running the simulations, you should run the initialization. During this phase the program will install a Python Virtual Environment with all the required packages
**Note**: This step should be done only one time after the installation
To run the initialization, go to the directory where you have downloaded the gwskysim code and run the init script:

.. code-block:: console

 cd gwskysim
 ./init-gwskysim.sh

This will create a directory called venv-gwskysim-py3 at the same level of the gwskysim directory.
It also will create a setup-gwskysim.sh setup script in the root directory of the package.

Setup
===================
Once the initialization is done, every time you want to run gwskysim scripts, you should run the setup script that you can find in the main root

.. code-block:: console

 source setup-gwskysim.sh

At the end of the script you will see a message
.. code-block:: console

 gwskysim setup done!

You are now ready to go. Enjoy!
