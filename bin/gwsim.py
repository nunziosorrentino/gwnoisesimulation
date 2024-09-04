#!/usr/bin/env python

"""

 Script name: gwskysim

 Description:
 Python script that creates GW simulated data, projected on existing or 
 future detectors, given the source parameters files, start GPS or date 
 time and the duration of the simulation. All the other optional arguments
 are listed in the help message.

"""
import json
import os
import platform
import sys
import time
from optparse import OptionParser
import pandas as pd
import numpy as np
import h5py
import multiprocessing

from gwdama.io import GwDataManager
# custom imports
from gwskysim.utilities.gwlogger import GWLogger
# simualtion imports
from gwskysim.sources.gwpopulation import GWPopulation
from gwskysim.sources.gwsignal import GWSignal

# Now, some additional information

__author__ = "Nunziato Sorrentino"
__copyright__ = "Copyright 2019-2020, Nunziato Sorrentino"
__credits__ = ["Line for credits"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "M. Razzano"
__developer__ = "N. Sorrentino"
__email__ = ["nunziato.sorrentino@pi.infn.it", "massimiliano.razzano@pi.infn.it"]
__status__ = "Production"

# General variables
os_system = platform.system()
running_python_version = int(sys.version[0])
work_dir = os.getcwd()
script_name = os.path.split(sys.argv[0])[1].split('.')[0]
script_path = os.path.join(os.getcwd(), os.path.dirname(sys.argv[0]))

# General Functions. None here, already imported from GeneralUtilities

#######################################################
# Main
#######################################################

if __name__ == '__main__':

    usg = "\033[1;31m%prog [ options] FT1FileName\033[1;m \n"

    desc = "\033[34mPrepare the files containing sources for the gw simulation\033[0m"

    parser = OptionParser(description=desc, usage=usg)
    parser.add_option("-s", "--starttime", type=float, help="Start GPS or date time of the")
    parser.add_option("-d", "--duration", type=float, help="Duration of the simulation in seconds")
    parser.add_option("-r", "--samplerate", type=int, help="Sample rate of the data set in Hz")
    parser.add_option("-i", "--inter", default=False, help="Interactive mode, showing strain plot") 
    parser.add_option("--detectors", default="V1 H1 L1", type=str,
                      help="List of detectors at which the wave projection is applicated, in the string spaced format:\n"+\
                           "if I want to project to LIGO/Virgo detector, the string is 'V1 H1 L1'")
    parser.add_option("--approximant", default=None, type=str,
                      help='Approximant used for the templates configurations')
    parser.add_option("--lowerfreq", default=10, type=float,
                      help='Lower frequency of the waveform in case of cbc objects')
    parser.add_option("--minsnr", default=1, type=float,
                      help='Minimum SNR value for the transient sources simulation')
    parser.add_option("--noise", choices=[None, 'white', 'colored'], default=None,
                      help="If not None, choice between 'white' and 'colored'") 
    parser.add_option("--nfiles", default=1, type=int,
                      help="Number of files, each grouping detectors ones, given in the simulator output")                                             
    parser.add_option("--nproc", default=None, type=int,
                      help='If needed, select the number of processes for using multiprocessing')                  
    parser.add_option("--fileformat", choices=['gwdama', 'gwosc'], default='gwdama',
                      help="Format of the output hdf file, choice between ['gwdama', 'gwosc']")                  
    parser.add_option("--outfileprefix", default=None, type=str,
                      help='Prefix of the output hdf file(s)')                                                                                        
    parser.add_option("--outdir", default=None, type=str, help="Output Directory")
    parser.add_option("--verbose", default=50, type=int, 
                      help="Verbose: every how many sources you want to be notified about the simulation status?")                 
    parser.add_option("--debug", default=False, help="Debug mode")

    (options, args) = parser.parse_args()

    start_time = options.starttime
    duration = options.duration
    sample_rate = options.samplerate
    detectors = options.detectors
    approximant = options.approximant
    lower_freq = options.lowerfreq
    min_snr = options.minsnr
    noise = options.noise
    n_files = options.nfiles
    n_proc = options.nproc
    file_format = options.fileformat
    if options.outfileprefix is None:
        of_prefix = 'gwskysim_{}'.format(file_format)
    else:
        of_prefix = options.outfileprefix    
    output_dir = options.outdir
    verbose = options.verbose
    if options.inter=='True':
        plot = True
    else:
        plot = False
    debug = options.debug

    detectors = detectors.split()
    population_files = ()

    #####################################################

    # Start the logger and check if debug mode is on
    my_logger = GWLogger("generate_pop_logger")
    my_logger.set_loglevel("INFO")

    if debug:
        my_logger.set_loglevel("DEBUG")

    my_logger.info('**************************************')
    my_logger.info("    " + str(script_name))
    my_logger.info("    (Running on " + os_system + " OS)")
    my_logger.info('**************************************')

    if len(args) == 0:
        my_logger.fatal("No HDF5 parameters file specified!")
        my_logger.fatal('Type ' + str(script_name) + '.py -h for help\n')
        sys.exit(1)
    elif len(args)==1 and args[0]=='None':
        my_logger.info("None argument means only background component is required for simulation!")    
    else:
        population_files = args
        
    # List parameters and args
    my_logger.info('\nInput options:')
    for key, val in iter(parser.values.__dict__.items()):
        my_logger.info(key + ": " + str(val))
        
    if output_dir is None:
        output_dir = ''
        for i in os.getenv('PYTHONPATH').split(':'):
            if i.endswith('gwskysim'):
                output_dir=os.path.abspath(os.path.join(i, '..', 'GWsimdata'))
                
    my_logger.debug("Output directory set to {}".format(output_dir))

    if not os.path.exists(output_dir):
        os.mkdir(output_dir) 

    
    src_population = GWPopulation(*population_files)
    my_logger.info(src_population) 
    n_sources = len(src_population)
    
    noise_dict = {'No Noise': None,
                  'White Noise': 'white',
                  'Colored Noise': 'colored'}

    new_duration = duration/n_files  
    new_start_times = np.linspace(start_time, start_time+duration, 
                                  n_files+1)[:-1]
                                  
    #strain_signals_list = []                              
    for i_f in range(n_files):
        if n_proc is None:
            strain_signals = src_population.generate_waveforms(new_start_times[i_f], 
                                                               new_duration, 
                                                               sample_rate, 
                                                               detectors=detectors,
                                                               approximant=approximant, 
                                                               lower_freq=lower_freq,
                                                               #min_snr=min_snr,
                                                               verbose=verbose,
                                                               noise=noise)
    
            # Save strain signals in hdf5 files
            my_logger.info("Saving hdf file(s): {}/{}".format(len(detectors)*(i_f+1), 
                                                        len(detectors)*n_files))
            for single_s in strain_signals:                 
                single_s.save(file_format, output_dir, of_prefix)
    
        else:       
            def generate_waveforms(iminmax):
                selected_pop = src_population[iminmax[0]:iminmax[1]]
                if iminmax[0] == 0:
                    noise_n = noise
                if iminmax[0] != 0:
                    noise_n = None
                signals_ = selected_pop.generate_waveforms(new_start_times[i_f],
                                                           new_duration, 
                                                           sample_rate, 
                                                           detectors=detectors,
                                                           approximant=approximant, 
                                                           lower_freq=lower_freq,
                                                           #min_snr=min_snr,
                                                           verbose=verbose,
                                                           noise=noise_n)
                return signals_
        
            # defining parallel processes for generating waveforms
            func_inputs = []
            for x in range(n_proc):
                imin = int(x*n_sources/n_proc)
                imax = int((x+1)*n_sources/n_proc)
                func_inputs.append([imin, imax])
            p = multiprocessing.Pool(processes=n_proc)
            proc_outputs = p.map(generate_waveforms, func_inputs)
            p.close()
            p.join()

            noise_string = list(noise_dict.keys())[list(noise_dict.values()).index(noise)]
            strain_signals = [GWSignal(noise_string, new_start_times[i_f], new_duration, sample_rate, 
                                       detector=d) for d in detectors]   
            #strain_signals[0].plot()
            for s_ in proc_outputs:
                for idet in range(len(detectors)):
                    strain_signals[idet] += s_[idet]
                    stri_list = [stri for stri in s_[idet].name.split() if stri[0].isdigit()]
                    strain_signals[idet].name += " "+" ".join(stri_list)
                    #strain_signals[0].plot()                                     
            
            # function for generating output hdf5 files
            def generate_hdf5(idet):
                strain_signals[idet].save(file_format, output_dir, of_prefix)

            my_logger.info("Saving hdf file(s): {}/{}".format(len(detectors)*(i_f+1), 
                                                        len(detectors)*n_files))
            # parallel processes for generating output hdf5 files
            p = multiprocessing.Pool(processes=len(detectors))
            list_idet = list(range(len(detectors)))
            p.map(generate_hdf5, list_idet)
            p.close()
            p.join()  
        #strain_signals_list.append(strain_signals)  
        del strain_signals
        where_saved = os.path.join(output_dir, '{}_{}_{}-{}GPS.hdf5'.format(of_prefix, '{DET}', 
                                   int(new_start_times[i_f]), int(new_start_times[i_f]+new_duration)))
        my_logger.info('All waveforms have been generated and saved into {}'.format(where_saved))        
                 
    if plot:
        import glob
        import matplotlib.pyplot as plt
        hdf_files = glob.glob(os.path.join(output_dir, '{}_*_*-*GPS.hdf5'.format(of_prefix)))
        for hfile in hdf_files:
            new_dama = GwDataManager(hfile, mode='r')
            new_dset = new_dama['gwskysim_simulation']
            new_dset.to_TimeSeries().plot()
            plt.title(new_dset.attrs['channel'])
            #plt.savefig('{}_{}.png'.format(of_prefix, new_dset.attrs['channel'][:2]))
            new_dama.close()
            del new_dama
        plt.show()
                             
   
