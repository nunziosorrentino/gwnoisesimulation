import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from gwskysim.sources.gwbinarysource import GWBinarySource
from gwskysim.sources.gwpopulation import GWPopulation
from gwskysim import GWSKYSIM_CONF
from gwdama.io import GwDataManager
from pycbc.filter import matched_filter

stype = 'BHBH-CE'
file_base = 'gws000000-gws000999'
duration = 1000
TIME = []
SNRs = [[], [], []]
dcts = ["V1", "L1", "H1"]

#pop_files_base = os.path.join(GWSKYSIM_CONF, '..', '..')
s_files = {}
for d in dcts:
    sign_file = os.path.join('..', '..',  'GWsimdata', '{}_gwskysim-simulation_{}_1242103936-1242104936GPS.hdf5'.format(stype, d))
    s_files[d] = sign_file
print("S. FILEs:", s_files)	
pop_files = os.path.join('..', '..',  'GWsim_{}_20220629'.format(stype), 'sources_parameters', '{}.hdf5'.format(file_base))
pop_df = pd.read_hdf(pop_files, "df_sources_parameters")
samples = len(pop_df)

plt.figure()
counts, _, _ = plt.hist(pop_df['distance'], bins=50)

a = 3
x_ = np.linspace(0, 1, 100)
y_ = (a+1)*x_**a
norm_y = y_*(counts.max()/y_.max())
norm_x = x_*(pop_df['distance'].max()-pop_df['distance'].min())+pop_df['distance'].min()
plt.plot(norm_x, norm_y, label='power{}'.format(a))
plt.xlabel("Distance [pc]")
plt.title("{}".format("BHBH-CE"))
plt.legend()
plt.show()

src_pop = GWPopulation(pop_files)

psd_dict = {}
for dct in dcts:
    asd_file = os.path.join(GWSKYSIM_CONF, "{}_O3.txt".format(dct))
    f_asd, asd = np.loadtxt(asd_file, unpack=True) 
    psd_dict[dct] = (f_asd, asd**2)

fs = 4096
for i in range(len(src_pop)):
    signal = src_pop[i].generate_signal(src_pop[i].tp-0.5,
                                        1, fs, 
                                        detectors=dcts)  
                                     
    #signal_noise = src_pop.generate_waveforms(src_pop[i].tp-0.5, 1, fs,
    #                                   detectors=dcts, noise='colored')                               
                                                                 
    TIME.append(src_pop[i].tp)                           
    for di in range(len(dcts)):
        
        with GwDataManager(s_files[dcts[di]], mode='r') as strain_dama: 
            strain_dset = strain_dama['gwskysim_simulation']
            f_, psd_ = strain_dset.psd(8, 4, return_type='array') 
            #plt.figure()
            #plt.loglog(f_.value, np.sqrt(psd_.value), label=dcts[di])
            #plt.legend()
            #plt.xlim((10, 1000))
            #plt.show()
            times_ = np.linspace(strain_dset.attrs["t0"], strain_dset.attrs["t0"]+duration, len(strain_dset))                            
            signal_noise = strain_dset[(times_>=src_pop[i].tp-0.5)&(times_<=src_pop[i].tp+0.5)]
        sampled_times = np.linspace(src_pop[i].tp-0.5, src_pop[i].tp+0.5, fs) 
        print("DETECTOR: ", dcts[di])
        snr_v = GWBinarySource.SNR(sampled_times, 
                                   signal_noise, 
                                   signal[di].get_sampled_time(), 
                                   signal[di].sampled_strain, 
                                   psd_, t0=src_pop[i].tp)
        print("SNR =", snr_v)
        #signal[di].plot()
        #signal_noise.plot()
        #plt.show()
        SNRs[di].append(snr_v)  
    print('Source number: {}'.format(i+1))                           

plt.figure()
counts_v, _, _ = plt.hist(SNRs[0], bins=10, label='V1')
counts_l, _, _ = plt.hist(SNRs[1], bins=10, label='L1')
counts_h, _, _ = plt.hist(SNRs[2], bins=10, label='H1')
plt.xlabel("SNR")
plt.title("{}".format("BHBH-CE"))
plt.legend()

plt.show()

mt_df = pd.DataFrame({"TIME":TIME, "SNR_v1":SNRs[0], "SNR_l1":SNRs[1], "SNR_h1":SNRs[2]})                         
mt_df.to_hdf(os.path.join('..', '..', 'GWsimdata', '{}_labels_{}.hdf5'.format(stype, file_base)), '{}'.format(stype))                     
print(mt_df)                                                  
