# -*- coding: utf-8 -*-
"""
Created on Sat Jun 30 11:24:23 2018

@author: will
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, lfilter

def cosine(t, amp, freq, phase = 0):
    return amp*np.cos(2*np.pi*freq*t + phase)

def view_section(datax, datay, samp_density, osc_win, Nosc = 10):
    #Just take a section for viewing: first Nosc oscillations
    Psamp = samp_density/osc_win #number of samples per period
    sec = int(Psamp*Nosc)
    return datax[0:sec], datay[0:sec]

def no_off_spectrum(sig, t):
    f = np.fft.rfftfreq(len(sig), d = t[1]-t[0])
    sig -= np.mean(sig)
    y = np.fft.rfft(sig)
    return f, y

class anti_alias:
    def __init__(self, tdata, ydata):
        self.ydata = ydata
        self.tdata = tdata

    def butter_bandpass(self, lowcut, highcut, fs, order=5):
        nyq = 0.5*fs
        low = lowcut/nyq
        high = highcut/nyq
        b, a = butter(order, [low, high], btype='band')
        return b, a

    def butter_bandpass_filter(self, lowcut, highcut, fs, order=5):
        b, a = self.butter_bandpass(lowcut, highcut, fs, order=order)
        y = lfilter(b, a, self.ydata)
        return y

    def sig_shift(self, freq):
        return self.tdata*freq

    def sample(self, samp_f):
        return self.ydata[0::samp_f]

    def aa_filt(self, freq, fs):
        yfilt = self.butter_bandpass_filter(freq - 0.1*freq, freq + 0.1*freq, fs)
        yshift = self.sig_shift(freq)
        #STILL NEED TO FINISH THIS!! MAKE SURE YOU DON'T MIX UP X's and Y's

def main():
    #define some parameters
    Amp = 2.
    freq = 50.
    nyq = 2*freq
    window = 100
    osc_win = window*freq #number of oscillations in the window

    samp_dense = 1e6
    t = np.linspace(0, window, samp_dense)
    y = cosine(t, Amp, freq)
    t_sec, y_sec = view_section(t, y, samp_dense, osc_win)

    #samp_sparse
    eps = 0.01*osc_win
    samp_sparse = 1.9*osc_win
    #samp_sparse = 2*osc_win - eps
    ts = np.linspace(0, window, samp_sparse)
    ys = cosine(ts, Amp, freq)
    ts_sec, ys_sec = view_section(ts, ys, samp_sparse, osc_win)

    #Now use anti-aliasing BPF before sampling
    T = 1/freq
    fs = samp_dense/T
    filt = anti_alias(t, y)



    #Take some FFTs
    fd, yd = no_off_spectrum(y, t)
    fs, ys = no_off_spectrum(ys, ts)

    #plot formatting
    plt.rc('text', usetex = False)
    plt.rc('font', family = 'serif')
    plt.rc('font', size = 22)
    plt.rc('axes', linewidth = 2)
    plt.rc('lines', linewidth = 3)
    plt.rc('legend', fontsize = 16)
    plt.rc('figure', figsize = (10, 6))
    plt.rc('lines', markersize = 15)

    #plots
    plt.figure(1)
    plt.plot(t_sec, y_sec)
    plt.plot(ts_sec, ys_sec, '.', color = 'red')
    plt.ylim(-2*Amp, 2*Amp)

    #plot FFTs
    plt.figure(2)
    plt.plot(fd, abs(yd)/max(abs(yd)))
    plt.plot((fs-1.9*freq)*-1, abs(ys)/max(abs(ys)))
    plt.xlim(-freq-eps, freq+eps)
    plt.ylim(-0.5, 1.5)

if __name__ == "__main__":
    main()