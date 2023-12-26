import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__),'../..'))
import numpy as np
from scipy.signal import stft


def GW_spectrogram(sim_dim, GWs, window_size):
    assert sim_dim in [1, 2, 3], "Simulation MUST be 1, 2 or 3D."
    if sim_dim == 1:
        return GW_spectrogram_1D(GWs)
    elif sim_dim == 2:
        return GW_spectrogram_2D(GWs, window_size)
    else:
        return GW_spectrogram_3D(GWs, window_size)

def GW_spectrogram_1D(GWs):
    print("And also no spectrogram for you :'(\nごめんなさい")
    return None

def GW_spectrogram_2D(GWs, window_size):
    """
        
    """
    fs = (GWs[1, 0] - GWs[0, 0]) ** (-1)
    frequency, time, Zxx = stft(GWs[:,1], fs = fs, nperseg=window_size)
    time = time / time[-1] * GWs[-1, 0]
    return time, frequency, np.abs(Zxx)

def GW_spectrogram_3D(GWs, window_size):
    fs = (GWs[1, 0] - GWs[0, 0]) ** (-1)
    freq_pl_e, time_pl_e, Zxx_pl_e = stft(GWs[:,1], fs = fs, nperseg=window_size)
    time_pl_e = time_pl_e / time_pl_e[-1] * GWs[-1, 0]
    freq_pl_p, time_pl_p, Zxx_pl_p = stft(GWs[:,2], fs = fs, nperseg=window_size)
    time_pl_p = time_pl_p / time_pl_p[-1] * GWs[-1, 0]
    freq_cr_e, time_cr_e, Zxx_cr_e = stft(GWs[:,3], fs = fs, nperseg=window_size)
    time_cr_e = time_cr_e / time_cr_e[-1] * GWs[-1, 0]
    freq_cr_p, time_cr_p, Zxx_cr_p = stft(GWs[:,4], fs = fs, nperseg=window_size)
    time_cr_p = time_cr_p / time_cr_p[-1] * GWs[-1, 0]
    
    return np.stack((time_pl_e, time_pl_p, time_cr_e, time_cr_p), axis=-1), \
           np.stack((freq_pl_e, freq_pl_p, freq_cr_e, freq_cr_p), axis=-1), \
           np.stack((np.abs(Zxx_pl_e), np.abs(Zxx_pl_p), np.abs(Zxx_cr_e),
                     np.abs(Zxx_cr_p)), axis=-1)