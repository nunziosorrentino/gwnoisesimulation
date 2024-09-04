import numpy as np
import os
import pandas as pd
import json
import logging
import matplotlib.pyplot as plt
from gwskysim.sources.gwsignal import GWSignal

from pycbc.waveform import get_td_waveform

from numpy.fft import irfft, rfftfreq
from scipy import integrate
from scipy.interpolate import interp1d
from pycbc.types.timeseries import TimeSeries
from pycbc.types.frequencyseries import FrequencySeries
from pycbc.filter import matched_filter

class GWSource:
    """Base class for all the model components (i.e., the GW sources).

        Given a certain time interval and sample frequency, the idea is to
        generate the equivalent of the stain amplitude in a single detector
        through an object incapsulating all the relevant information
        (e.g. source parameters, time of emission (if transient), average
        noise spectrum)

        Arguments
        ---------
        name : string
        The name of the source.
    """

    def __init__(self, name):
        """Constructor.
        """
        self.name = name

    @property
    def parameters_dict(self):
        """Property method that generate a dictionary with source
           parameters.
        """
        return {i: self.__dict__[i] for i in self.__dict__.keys() \
                if i != 'name'}

    @property
    def parameters_df(self):
        """Property method that generate a managefull source parameters
           data frame.
        """
        return pd.DataFrame(self.parameters_dict, index=[0])

    @classmethod
    def from_file(cls, m_source_file):
        """Class method that reads parameters from a signal file.

        Warning
        -------
        Usable only with JSON files, other types will be implemented
        in next versions.
        """
        # Assert if the file exist and if it is a JSON file
        if not os.path.isfile(m_source_file):
            raise FileNotFoundError("Please verify your input parameters file!!")

        if m_source_file.split('.')[1] == 'json':
            # Read the JSON file
            with open(m_source_file) as json_file:
                kys_vls = json.load(json_file)

        if m_source_file.split('.')[1] == 'csv':
            kys_vls = pd.read_csv(m_source_file)

        if m_source_file.split('.')[1] == 'hdf5':

            kys_vls = pd.read_hdf(m_source_file)

        else:
            raise NotImplementedError("Only JSON files can be uploaded.")

        return cls(**kys_vls)

    @classmethod
    def from_df(cls, data_frame_r):
        """Class method that reads parameters from a pandas dataframe raw.

           Parameters
           ----------
           data_frame_r : `~pandas.core.series.Series` instance
           Raw of the dataframe containing all the parameters defining a
           single source.

        """
        kys_vls = dict(data_frame_r.dropna())
        kys_vls['name'] = kys_vls.pop('id')

        del kys_vls['type']

        return cls(**kys_vls)

    @staticmethod
    def SNR(t_stretch, stretch, t_template, template, psd, t0=0):
           # Get fft of data stretch
           dt = t_stretch[1]-t_stretch[0]
           stretch_ts = TimeSeries(stretch, dt, t_stretch[0])
           stretch_fs = stretch_ts.to_frequencyseries(0.33)#0.83)
           #print('AAA', stretch_fs.get_delta_f())
           # Get fft of template
           dt_i = abs(t_stretch[0] - t_template[0])
           dt_f = abs(t_stretch[-1] - t_template[-1])
           n_i = int(dt_i*(1/dt))+1
           n_f = int(dt_f*(1/dt))+1
           template_ = np.pad(template, (n_i, n_f), "constant", constant_values=(0, 0))
           template_ts = TimeSeries(template_, dt, t_stretch[0])
           template_fs = template_ts.to_frequencyseries(delta_f=stretch_fs.get_delta_f())
           # Get PSD
           #plt.figure()                     
           #plt.plot(template_fs.sample_frequencies, abs(template_fs))
           #plt.show()
           psd = FrequencySeries(psd, stretch_fs.get_delta_f(), t_stretch[0])
           #print(len(template_fs), len(stretch_fs))
           SNR_r = matched_filter(template_fs, template_fs, psd=psd,
                                  low_frequency_cutoff=10)
           #print("qui", SNR_r)   
           #ntemplate_ts = TimeSeries(template_/np.sqrt(abs(SNR_r))._data, dt, t_stretch[0])
           #ntemplate_fs = template_ts.to_frequencyseries(delta_f=stretch_fs.get_delta_f())                    
           SNR = matched_filter(template_fs, stretch_fs, psd=psd,
                                low_frequency_cutoff=10)
           t_snr = SNR.get_sample_times()
           snr_mask = [(t_snr>=t0-0.2)&(t_snr<=t0+0.2)]
           SNR = SNR._data[snr_mask]                   
           #plt.figure()                     
           #plt.plot(abs(SNR))
           #plt.show()
           return max(abs(SNR))      

    def _import_template(self, times_array, approximant=None, lower_freq=None,
                               sample_rate=None):
        r"""Return waveforms amplitude arrays 

            (i.e. $h_{+}$ and $h_{\times}$).

        Arguments
        ---------
        times_array : array of floats
            the times array of which one wants to evaluate hp and hc.


        approximant : str
            keyword argument required for the model import from used by
            pycbc bank.
        lower_freq : float
            Lower frequency required for the configuration of the 
            pycbc model in Hz.
        sample_rate : float
            Sampe rate value in Hz, just needed if data aren't equally time separated.            


        Note
        ----
        This use a __call__ method that should be re-implemented in every
        derived class.
        """

        if len(times_array) < 2:
            raise ValueError('Time Series must contain at least two samples.')
            
        if sample_rate is None:
            sample_rate = round(1/(times_array[1]-times_array[0]))            

        if approximant is None or approximant.startswith('GW'):
            if self(times_array) is None:
                return (None, None)
            try:
                hp_a, hc_a, time0 = self(times_array)
            except ValueError:
                hp_a, hc_a = self(times_array)
                time0=times_array[0]
            
            d_ = len(hp_a)/sample_rate  
                
            hp = GWSignal(self.name, time0, d_, sample_rate,
                          sampled_strain=hp_a, detector='')
            hc = GWSignal(self.name, time0, d_, sample_rate,
                          sampled_strain=hc_a, detector='')

        else:
            assert (type(approximant) == str)
            pycbc_dict = {i: self.parameters_dict[i] \
                          for i in self.parameters_dict.keys() \
                          if i!='ra' if i!='dec'\
                          if i!='polarization' if i!='distance'}             
            hp_s, hc_s = get_td_waveform(approximant=approximant, 
                                     delta_t=1/sample_rate,
                                     f_lower=lower_freq,
                                     **pycbc_dict)                 
            if 'tmerge' in self.parameters_dict.keys():
                merge_t = self.parameters_dict['tmerge']
            else:
                try:
                    np.random.seed(int(self.name[-6:]))
                except ValueError:
                    np.random.seed(abs(int(times_array.min())))    
                merge_t = np.random.randint(int(times_array.min()),
                                            int(times_array.max()))
                merge_t += np.random.uniform(0, 1)
                #merge_t = times_array.max()
            hp_s.start_time += merge_t            
            hc_s.start_time += merge_t
            
            hp = GWSignal(self.name, hp_s.start_time, hp_s.duration, 
                          hp_s.sample_rate, sampled_strain=hp_s, detector='')
            hc = GWSignal(self.name, hc_s.start_time, hc_s.duration,
                          hc_s.sample_rate, sampled_strain=hc_s, detector='')

        return hp, hc

    def generate_signal(self, m_start_time, m_duration, m_sample_rate,
                        **kwargs):
        """Generate the strain signal in the proper detector.

        Note
        ----
        This is a do-nothing function that will be implemented soon.
        """
        raise NotImplementedError

    def __call__(self, times_array):
        """This is a __call__ method that should be re-implemented in every
           sub-class of GWSource.
        """
        raise NotImplementedError

    def __str__(self):
        """String formatting.
        """
        return '{} "{}": \n{}'.format(self.__class__.__name__, self.name,
                                      self.parameters_df.iloc[0])

"""
    @staticmethod
    def gluer(signal_amplitude, signal_times, signal_t0, total_signal_times):
        ##incolla una serie temporale signal times,centrata all'istante t0 rispetto alla serie temporale
        # totale total_signal_times
        # dt = np.abs(signal_times[1] - signal_times[2])
        total_samples = len(total_signal_times)
        total_duration = total_signal_times[total_samples - 1] - total_signal_times[0]
        index = np.int((signal_t0 * total_samples) / total_duration)
        # print("index",index)
        window_lenght = len(signal_times)  # lenght of the signal
        window_starting_index = np.int(index - window_lenght / 2)
        # print("init_index",window_starting_index)
        window_end_index = window_starting_index + window_lenght
        aux = np.zeros(total_samples)

        aux1 = np.append(aux[0:window_starting_index], signal_amplitude)

        aux2 = np.append(aux1, aux[window_end_index:])
        # print((signal_amplitude))
        return aux2
"""        

