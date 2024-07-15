#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 20:20:46 2024

@author: yaozhiyan
"""

import sys
import h5py
import MeasurementInfo
from MeasurementInfo import MeasurementInfo
from RunInfo import RunInfo
from scipy import signal
import ProcessWaveforms_MultiGaussian
from ProcessWaveforms_MultiGaussian import WaveformProcessor as WaveformProcessor
from AnalyzePDE import SPE_data
from AnalyzePDE import Alpha_data
import numpy as np
import matplotlib.pyplot as plt
from array import array
import pickle

files = ['Run_1719510386.hdf5', 'Run_1719511233.hdf5', 'Run_1719511764.hdf5', 'Run_1719512545.hdf5', 'Run_1719512918.hdf5',
         'Run_1719513299.hdf5', 'Run_1719513664.hdf5', 'Run_1719514159.hdf5', 'Run_1719514436.hdf5']

upperlim = [0.06,
            0.15,
            0.05,
            0.14,
            0.14,
            0.29,
            0.06,
            0.14,
            0.05]

proms = [0.005, 0.005, 0.005, 0.005, 0.006,
         0.008, 0.005, 0.006, 0.005]



# %%
run_spe_solicited = RunInfo(['Run_1719515422.hdf5'],
                            specifyAcquisition=False,
                            do_filter=True,
                            is_solicit=True,
                            upper_limit=4.4,
                            baseline_correct=True,
                            prominence=0.0001,
                            plot_waveforms=False)

# %%
runs = []
for file_ in range(len(files)):
    # print('This is proms', proms[file_])
    run_spe = RunInfo([files[file_]],
                      specifyAcquisition=False,
                      do_filter=True,
                      upper_limit=upperlim[file_],
                      baseline_correct=True,
                      prominence=proms[file_],
                      plot_waveforms=False)
    runs.append(run_spe)

# # %%
# with open('runs.pkl', 'wb') as f:
#     pickle.dump(runs, f)
# # #%%
# # print(file_)
# # print(len(runs))
# %%
for run in runs:
    run.plot_hists("0", "0", new=True)
# %%
biases = [run.bias for run in runs]
cutoffs = []
centers_guesses = []
for run in runs:
    bins = int(round(np.sqrt(len(run.all_peak_data))))
    # print('THis is len',len(run.all_peak_data))
    count, edges = np.histogram(run.all_peak_data, bins=bins)
    print('THis is bias', run.bias)
    if run.bias == 33:
        d = 5
        prom = 10
    if run.bias == 33.5:
        d = 5
        prom = 15
    if run.bias == 34.5:
        d = 5
        prom = 25
    if run.bias == 35.5:
        d = 2
        prom = 15
    if run.bias == 35:
        d = 2
        prom = 15
    if run.bias == 36.5:
        d = 2
        prom = 10
    if run.bias == 37:
        d = 5
        prom = 20
    if run.bias == 36:
        d = 2
        prom = 10
    if run.bias == 32:
        d = 10
        prom=20
    if run.bias == 32.5:
        d = 10
        prom=30
    if run.bias == 34:
        d = 5
        prom=25
    # else:
    #     print(runs.index(run.bias))

        
    #     d = 10
    #     prom = 25


    centers = (edges[:-1]+edges[1:])/2
    peaks, props = signal.find_peaks(count, prominence=prom, distance=d)

    # print('This is peaks', peaks)

    fitrange = ((centers[peaks[3]]-centers[peaks[0]])/2)

    if run.bias == 34.5:

        range_low = centers[peaks[0]]-0.2*fitrange
        range_high = centers[peaks[3]]+0.55*fitrange
    if run.bias == 37:

        range_low = centers[peaks[0]]-0.3*fitrange
        range_high = centers[peaks[3]]+0.75*fitrange
    if run.bias == 35.5:

        range_low = centers[peaks[0]]-0.2*fitrange
        range_high = centers[peaks[3]]+0.5*fitrange
    if run.bias == 36.5:

        range_low = centers[peaks[0]]-0.3*fitrange
        range_high = centers[peaks[3]]+0.75*fitrange
    if run.bias == 33.5:

        range_low = centers[peaks[0]]-0.3*fitrange
        range_high = centers[peaks[3]]+0.75*fitrange
    if run.bias == 36:

        range_low = centers[peaks[0]]-0.3*fitrange
        range_high = centers[peaks[3]]+0.9*fitrange
    if run.bias == 35:

        range_low = centers[peaks[0]]-0.3*fitrange
        range_high = centers[peaks[3]]+0.6*fitrange
    if run.bias == 32:

        range_low = centers[peaks[0]]-0.4*fitrange
        range_high = centers[peaks[3]]+0.8*fitrange
    if run.bias ==33:
        range_low = centers[peaks[0]]-0.3*fitrange
        range_high = centers[peaks[3]]+0.6*fitrange
    if run.bias == 32.5:

        range_low = centers[peaks[0]]-0.3*fitrange
        range_high = centers[peaks[3]]+0.6*fitrange
    if run.bias ==34:
        range_low = centers[peaks[0]]-0.3*fitrange
        range_high = centers[peaks[3]]+0.5*fitrange
    # else:

    #     range_low = centers[peaks[0]]-0.6*fitrange
    #     range_high = centers[peaks[3]]+0.3*fitrange
        
    peaks = peaks[0:]
    cutoffs.append((range_low, range_high))

    centers_guesses.append([centers[peak] for peak in peaks])

    plt.yscale('log')
    plt.hist(np.array(run.all_peak_data), bins=bins)
    plt.plot(centers[peaks], count[peaks], 'r.')

    plt.show()

# %%
for n in range(len(runs)):
    T = 169.5
    con = 'LXe'
    info = MeasurementInfo()
    info.condition = con
    info.date = runs[n].date
    info.temperature = T
    info.bias = biases[n]
    info.baseline_numbins = 150
    info.peaks_numbins = 200
    info.data_type = 'h5'
    wp = WaveformProcessor(info=info,run_info_self=runs[n],run_info_solicit=run_spe_solicited, centers=centers_guesses[n], cutoff=cutoffs[n])
            
    wp.process(do_spe=True, do_alpha=False)
    wp.plot_peak_histograms(log_scale=False)
    wp.plot_spe()
# %%
invC = 0.0132
invC_err = 0.000008
# %% make a list of ProcessWaveforms objects
campaign = []
for i in range(len(runs)):
    info_spe = MeasurementInfo()
    info_spe.condition = 'LXe'
    info_spe.temperature = 170
    info_spe.bias = runs[i].bias
    info_spe.baseline_numbins = 50
    info_spe.peaks_numbins = 150
    info_spe.data_type = 'h5'
    wp = WaveformProcessor(info_spe, centers=centers_guesses[i],
                           run_info_self=runs[i],
                           run_info_solicit=run_spe_solicited,
                           baseline_correct=True, cutoff=cutoffs[i])
    wp.process(do_spe=True, do_alpha=False)
    campaign.append(wp)

spe = SPE_data(campaign, invC, invC_err, filtered=False)
# %%
spe.plot_spe(in_ov=False, absolute=False,
             out_file='/Users/yaozhiyan/Documents/Research/nEXO/20240627_SPE_LED_405_LXe/LXe.csv')

# %%
spe_CA = SPE_data(campaign, invC, invC_err, filtered=True)
spe_CA.plot_CA(
    out_file='/Users/yaozhiyan/Documents/Research/nEXO/20240627_SPE_LED_405_LXe/LXe_CA.csv')