import numpy as np
import os
import pylab
import h5py
import time
from scipy import interpolate
#from pycbc.types.timeseries import TimeSeries
from gwdama.io import GwDataManager

class GWSignal:
    def __init__(self, name, start_time, duration, sample_rate, 
                 sampled_strain=None, detector=None):
        """"Signal constructor."""
        if sampled_strain is None:
            sampled_strain = np.zeros(int(duration*sample_rate), dtype=int)
        #TimeSeries.__init__(self, initial_array=sampled_strain, 
        #                    delta_t=1/sample_rate, epoch=start_time)
        self.name = name                    
        self.start_time = float(start_time)
        self.duration = float(duration)
        self.end_time = self.start_time + self.duration                   
        self.sample_rate = int(sample_rate)
        self.delta_t = 1/self.sample_rate                              
        self.sampled_strain = np.array(sampled_strain)
        self.detector = detector
        
    def __getitem__(self, index):
    
        if isinstance(index, slice):
            s_sampled_strain = self.sampled_strain[index]
            
            s_start_time = self.start_time
            if index.start is not None:
                s_start_time += index.start*self.delta_t
            
            s_end_time = self.end_time    
            if index.stop is not None:
                s_end_time -= index.stop*self.delta_t
                
            s_duration = s_end_time - s_start_time
            return self.__class__(self.name, s_start_time, s_duration,
                                  self.sample_rate, s_sampled_strain,
                                  self.detector)
        else:
            return self.sampled_strain[index]        

    def __add__(self, other):
    
        assert(self.sample_rate==other.sample_rate)
        assert(self.detector==other.detector) 
        
        padd_min = int(abs(self.start_time-other.start_time)*self.sample_rate) 
        padd_max = int(abs(self.end_time-other.end_time)*self.sample_rate)

        sampled_strain_1 = self.sampled_strain
        
        # Create zero padding appropriately
        if self.start_time <= other.start_time and self.end_time >= other.end_time:
            #print("start time",self.start_time)
            #print("start time 222",other.start_time)
            sampled_strain_2 = np.pad(other.sampled_strain, 
                                     (padd_min, padd_max), mode='constant')
            # The rest of zeros to be added in case of wrong integer  
            # approximation in padd min and max calculation                          
            z_rest = len(sampled_strain_1) - len(sampled_strain_2)
            rest_pmin = z_rest // 2
            rest_pmax = z_rest - rest_pmin
            sampled_strain_2 = np.pad(sampled_strain_2, 
                                   (rest_pmin, rest_pmax), mode='constant')  
            #print('1 PASSED!')                                              
            
        if self.start_time > other.start_time and self.end_time < other.end_time:

            s = padd_min
            e = len(other) - padd_max
            sampled_strain_2 = other[s:e].sampled_strain
            # The rest of zeros to be added in case of wrong integer  
            # approximation in padd min and max calculation 
            z_rest = len(sampled_strain_1) - len(sampled_strain_2)
            rest_pmin = z_rest // 2
            rest_pmax = z_rest - rest_pmin
            r_s = rest_pmin
            r_e = len(other) - rest_pmax
            sampled_strain_2 = other[r_s:r_e].sampled_strain  
            #print('2 PASSED!')   
        
        if self.start_time < other.start_time and self.end_time < other.end_time:
            e = len(other) - padd_max
            sampled_strain_2 = other[:e].sampled_strain 
            sampled_strain_2 = np.pad(sampled_strain_2, 
                                     (padd_min, 0), mode='constant')
            if len(sampled_strain_2)>len(sampled_strain_1):
                diff_rim = len(sampled_strain_2)-len(sampled_strain_1)
                sampled_strain_2 = sampled_strain_2[diff_rim:]
            #print('3 PASSED!')                                                
                                     
        if self.start_time > other.start_time and self.end_time > other.end_time:
            s = padd_min
            sampled_strain_2 = other[s:].sampled_strain
            sampled_strain_2 = np.pad(sampled_strain_2, 
                                     (0, padd_max), mode='constant') 
            if len(sampled_strain_2)>len(sampled_strain_1):
                diff_rim = len(sampled_strain_2)-len(sampled_strain_1)
                sampled_strain_2 = sampled_strain_2[:-diff_rim]                         
            #print('4 PASSED!')                                                                                                                                     
        
        # Create the strain for the new signal
        new_sample_strain = sampled_strain_1 + sampled_strain_2                        
                     
        return self.__class__(self.name, self.start_time, 
                              self.duration, self.sample_rate, 
                              new_sample_strain, self.detector)  
                              
    def __iadd__(self, other):
        return self + other          
        
    def __len__(self):
        return len(self.sampled_strain)
        
    def __mul__(self, other):
  
        self.sampled_strain *= other
            
        return self     

    def get_sampled_time(self):
        """
        Get times array of the time series.
        """
        start_time = float(self.start_time)
        end_time = float(self.end_time)
        sampled_time = np.linspace(start_time, end_time,
                                   num=len(self))
        return sampled_time

    def to_gwdama(self, save=True, storage='tmp', 
                        m_output_prefix='gwskysim_dama'):
        """
        Method to write the dataset into a gwdama object, with the option
        of saving it in a hdf5 file gwdama formatted.
        
        Parameters
        ----------
        save : bool, True
            Save or not in a hdf5 file gwdama formatted;
        storage : str, tmp 
            Random data store ('mem') or temporaly in disk ('tmp');    
        m_output_prefix : str, gwskysim_dama
            Name of the output file prefix. Default output file name: 
            'gwskysim_dama_{detector}_{GPS Start}-{GPS Stop}GPS.hdf5';
        """
        m_dsetname = '{}_{}_{}-{}GPS'.format(m_output_prefix, self.detector, 
                                      int(self.start_time), int(self.end_time))  
                                                                 
        sim_dama = GwDataManager(m_dsetname, storage=storage)                                                        

        sim_dama.attrs["description"] = 'This is the result of a simulation made with gwskysim and saved in hdf5 format with gwdama!'
        #print('Inizio creazione dset')                                                                         
        sim_dset = sim_dama.create_dataset('gwskysim_simulation', 
                                           data=self.sampled_strain) 
        #print('Fine creazione dset')                                                            
        sim_dset.attrs['t0'] = float(self.start_time)
        sim_dset.attrs['unit'] = ' '
        sim_dset.attrs['channel'] = '{}:gwskysim_strain'.format(self.detector)
        sim_dset.attrs['sample_rate'] = int(self.sample_rate)
        sim_dset.attrs['name'] = self.name   
        
        if save:    
            sim_dama.write_gwdama(m_dsetname+'.hdf5')
            sim_dama.close()
            del sim_dama 
        if not save:
            return sim_dama                             
        
    def save_to_hdf(self, m_output_prefix='output_gwskysim'):
        """
        Method to write the dataset into an hdf5 file.
        
        Parameters
        ----------
        m_output_prefix : str, output_gwskysim
            Name of the output file prefix. Default output file name: 
            'output_gwskysim_{detector}_{GPS Start}-{GPS Stop}GPS.hdf5';
        """
        m_output_filename = '{}_{}_{}-{}GPS.hdf5'.format(m_output_prefix, 
                             self.detector, int(self.start_time), int(self.end_time))                                   

        if os.path.isfile(m_output_filename):
            raise FileExistsError('File already exists, no overwriting!')
        m_creation_time = str(time.strftime("%y-%m-%d_%Hh%Mm%Ss", 
                                                         time.localtime()))

        with h5py.File(m_output_filename, 'w') as m_out_hf:
            m_out_hf.attrs.create("time_stamp", m_creation_time)
            d_set = m_out_hf.create_dataset('gwskysim_simulation',
                                            data=self.sampled_strain, compression="gzip")
            d_set.attrs.create('GPS Start', float(self.start_time))
            d_set.attrs.create('Yunit', '')
            d_set.attrs.create('Xunit', 's')
            d_set.attrs.create('Channel', 'GWSignal')
            d_set.attrs.create('Sample Rate', int(self.sample_rate))
            d_set.attrs.create('Name', self.name)
            
    def save(self, f_format, m_out_dir, prefix=None):
        """"
        Save the signal strain in one of the format provided:
         - 'gwdama' like;
         - 'gwosc' like;
          
        Parameters
        ----------
        f_format : str
            specify here the format you prefer;
        m_out_dir : str
            specify here directory destination of the file;    
        prefix : str, None
            specify the prefix for your output file (if any);       
        """
        if f_format=='gwdama':
            if prefix is None:
                out_f = os.path.join(m_out_dir, 'gwskysim_dama')
                self.to_gwdama(save=True, m_output_prefix=out_f)
            else:
                out_f = os.path.join(m_out_dir, prefix)
                self.to_gwdama(save=True, m_output_prefix=out_f)
        elif f_format=='gwosc':
            if prefix is None:
                out_f = os.path.join(m_out_dir, 'output_gwskysim')
                self.save_to_hdf(m_output_prefix=out_f)
            else:
                out_f = os.path.join(m_out_dir, prefix)
                self.save_to_hdf(m_output_prefix=out_f) 
        else:
            raise AssertionError("Format not available, 'gwdama' and 'gwosc' like are available!")                      

    def plot(self):
        """
        Plot the signal strain time series.
        """
        pylab.figure(figsize=(20, 4))
        pylab.plot(self.get_sampled_time(), self.sampled_strain, 
                                      label=self.detector)
        pylab.title(self.name)
        pylab.ylabel('Strain []')
        pylab.xlabel('Time [s]')
        pylab.legend()
        #pylab.savefig(self.name)
        #pylab.close()

    def get_continue_strain(self):
        """
        """
        continue_strain = interpolate.interp1d(self.get_sampled_time(), 
                                      self.sampled_strain, kind='cubic')
        return continue_strain

    def smooth_start(self):
        """
        """
        len2smooth = int(self.sample_rate/16)
        hann_function = np.hanning(2*len2smooth)
        half_hann = hann_function[:len2smooth]
        smoothed_window = self.sampled_strain[:len2smooth] * half_hann
        self.sampled_strain[:len2smooth] = smoothed_window    
        
    def __str__(self):
        return "GWSignal '{}:{}' \n{}".format(self.detector, self.name, 
                                       self.sampled_strain.__str__())        

















