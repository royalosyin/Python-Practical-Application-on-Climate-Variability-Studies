"""
The following codes are extracted and modified from
Author:       Damien Irving, irving.damien@gmail.com
Description:  Calculate Fourier transform
"""

import math
import numpy as np
import scipy as sp
import scipy.fftpack as fftpack
from copy import deepcopy

def fourier_transform(signal, spacing):
    """Calculate the Fourier Transform.
    
    Args:
      signal (numpy.ndarray): Data to be transformed 
      spacing (scaler): sampling resolution
    
    Returns:
      sig_fft (numpy.ndarray): Coefficients obtained from the Fourier Transform
      freqs (numpy.ndarray): Wave frequency associated with each coefficient
    
    """        
    sig_fft = fftpack.fft(signal)
    sample_freq = fftpack.fftfreq(len(signal), d=spacing) * len(signal) * spacing  #units = cycles per length of domain
    sample_freq = np.resize(sample_freq, sig_fft.shape)
    
    return sig_fft, sample_freq
    
def spectrum(signal_fft, freqs, scaling='amplitude', variance=None):
    """Calculate the spectral density for a given Fourier Transform.
    
    Args:
      signal_fft, freqs (numpy.ndarray): Typically the output of fourier_transform()
      scaling (str, optional): Choices for the amplitude scaling for each frequency
        are as follows (see Wilks 2011, p440):
         'amplitude': no scaling at all (C)
         'power': sqaure the amplitude (C^2)
         'R2': variance explained = [(n/2)*C^2] / (n-1)*variance^2, 
         where n and variance are the length and variance of the 
         orignal data series (R2 = the proportion of the variance 
         explained by each harmonic)    
    """

    assert scaling in ['amplitude', 'power', 'R2']
    if scaling == 'R2':
        assert variance, \
        "To calculate variance explained must provide variance value" 
        
    if len(signal_fft.shape) > 1:
        print("WARNING: Ensure that frequency is the final axis")
    
    # Calculate the entire amplitude spectrum
    n = signal_fft.shape[-1]
    amp = np.abs(signal_fft) / n
    
    # The positive and negative half are identical, so just keep positive
    # and double its amplitude
    freq_limit_index = int(math.floor(n / 2)) 
    pos_amp = 2 * np.take(amp, range(1, freq_limit_index), axis=-1)
    pos_freqs = np.take(freqs, range(1, freq_limit_index), axis=-1)
    
    if scaling == 'amplitude':
        result = pos_amp
    elif scaling == 'power':
        result = (pos_amp)**2
    elif scaling == 'R2':
        result = ((n / 2) * (pos_amp**2)) / ((n - 1) * (variance))
    
    return result, pos_freqs
    
    
def inverse_fourier_transform(coefficients, sample_freq, 
                              min_freq=None, max_freq=None, exclude='negative'):
    """Inverse Fourier Transform.
    
    Args:
      coefficients (numpy.ndarray): Coefficients obtained from the Fourier Transform
      sample_freq (numpy.ndarray): Wave frequency associated with each coefficient
      max_freq, min_freq (float, optional): Exclude values outside [min_freq, max_freq]
        frequency range. (Note that this filtering keeps both the positive and 
        negative half of the spectrum)
      exclude (str, optional): Exclude either the 'positive' or 'negative' 
        half of the Fourier spectrum. (A Hilbert transform, for example, excludes 
        the negative part of the spectrum)
                                 
    """
    
    assert exclude in ['positive', 'negative', None]
    
    coefs = deepcopy(coefficients)  # Deep copy to prevent side effects
                                    # (shallow copy not sufficient for complex
                                    # things like numpy arrays)
    
    if exclude == 'positive':
        coefs[sample_freq > 0] = 0
    elif exclude == 'negative':
        coefs[sample_freq < 0] = 0
    
    if (max_freq == min_freq) and max_freq:
        coefs[np.abs(sample_freq) != max_freq] = 0
    
    if max_freq:
        coefs[np.abs(sample_freq) > max_freq] = 0
    
    if min_freq:
        coefs[np.abs(sample_freq) < min_freq] = 0
    
    result = fftpack.ifft(coefs)
    
    return result    
