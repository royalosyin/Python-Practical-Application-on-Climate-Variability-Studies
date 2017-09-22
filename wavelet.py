def wavelet(Y,dt,pad=0.,dj=0.25,s0=-1,J1=-1,mother="MORLET",param=-1):
    """
This function is the translation of wavelet.m by Torrence and Compo

import wave_bases from wave_bases.py

The following is the original comment in wavelet.m

#WAVELET  1D Wavelet transform with optional singificance testing
%
%   [WAVE,PERIOD,SCALE,COI] = wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
%
%   Computes the wavelet transform of the vector Y (length N),
%   with sampling rate DT.
%
%   By default, the Morlet wavelet (k0=6) is used.
%   The wavelet basis is normalized to have total energy=1 at all scales.
%
%
% INPUTS:
%
%    Y = the time series of length N.
%    DT = amount of time between each Y value, i.e. the sampling time.
%
% OUTPUTS:
%
%    WAVE is the WAVELET transform of Y. This is a complex array
%    of dimensions (N,J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
%    ATAN(IMAGINARY(WAVE),FLOAT(WAVE) gives the WAVELET phase.
%    The WAVELET power spectrum is ABS(WAVE)^2.
%    Its units are sigma^2 (the time series variance).
%
%
% OPTIONAL INPUTS:
% 
% *** Note *** setting any of the following to -1 will cause the default
%               value to be used.
%
%    PAD = if set to 1 (default is 0), pad time series with enough zeroes to get
%         N up to the next higher power of 2. This prevents wraparound
%         from the end of the time series to the beginning, and also
%         speeds up the FFT's used to do the wavelet transform.
%         This will not eliminate all edge effects (see COI below).
%
%    DJ = the spacing between discrete scales. Default is 0.25.
%         A smaller # will give better scale resolution, but be slower to plot.
%
%    S0 = the smallest scale of the wavelet.  Default is 2*DT.
%
%    J1 = the # of scales minus one. Scales range from S0 up to S0*2^(J1*DJ),
%        to give a total of (J1+1) scales. Default is J1 = (LOG2(N DT/S0))/DJ.
%
%    MOTHER = the mother wavelet function.
%             The choices are 'MORLET', 'PAUL', or 'DOG'
%
%    PARAM = the mother wavelet parameter.
%            For 'MORLET' this is k0 (wavenumber), default is 6.
%            For 'PAUL' this is m (order), default is 4.
%            For 'DOG' this is m (m-th derivative), default is 2.
%
%
% OPTIONAL OUTPUTS:
%
%    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
%           to the SCALEs.
%
%    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J1
%            where J1+1 is the total # of scales.
%
%    COI = if specified, then return the Cone-of-Influence, which is a vector
%        of N points that contains the maximum period of useful information
%        at that particular time.
%        Periods greater than this are subject to edge effects.
%        This can be used to plot COI lines on a contour plot by doing:
%
%              contour(time,log(period),log(power))
%              plot(time,log(coi),'k')
%
%----------------------------------------------------------------------------
%   Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
%
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made. This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%
% Notice: Please acknowledge the use of the above software in any publications:
%    ``Wavelet software was provided by C. Torrence and G. Compo,
%      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
%
% Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
%            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
%
% Please send a copy of such publications to either C. Torrence or G. Compo:
%  Dr. Christopher Torrence               Dr. Gilbert P. Compo
%  Research Systems, Inc.                 Climate Diagnostics Center
%  4990 Pearl East Circle                 325 Broadway R/CDC1
%  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
%  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
%----------------------------------------------------------------------------"""  
    #modules
    import numpy as np
    from wave_bases import wave_bases
    
    #set default
    n1 = len(Y)
    if (s0 == -1): s0=2.*dt
    if (dj == -1): dj = 1./4.
    if (J1 == -1): J1=np.fix((np.log(n1*dt/s0)/np.log(2))/dj)
    if (mother == -1): mother = 'MORLET'
    #print "s0=",s0
    #print "J1=",J1

    #....construct time series to analyze, pad if necessary
    x = Y - np.mean(Y);
    if (pad == 1):
        base2 = np.fix(np.log(n1)/np.log(2) + 0.4999)   # power of 2 nearest to N
        temp=np.zeros((2**(int(base2)+1)-n1,))
        x=np.concatenate((x,temp))
    
    n = len(x)

    #....construct wavenumber array used in transform [Eqn(5)]
    k = np.arange(1,np.fix(n/2)+1)
    k = k*(2.*np.pi)/(n*dt)
    k = np.concatenate((np.zeros((1,)),k, -k[-2::-1]));

    #....compute FFT of the (padded) time series
    f = np.fft.fft(x)    # [Eqn(3)]
    
    #....construct SCALE array & empty PERIOD & WAVE arrays
    scale=np.array([s0*2**(i*dj) for i in range(0,int(J1)+1)])
    period = scale.copy()
    wave = np.zeros((int(J1)+1,n),dtype=np.complex)  # define the wavelet array  # make it complex
    # loop through all scales and compute transform
    for a1 in range(0,int(J1)+1):
        daughter,fourier_factor,coi,dofmin=wave_bases(mother,k,scale[a1],param)
        wave[a1,:] = np.fft.ifft(f*daughter)  # wavelet transform[Eqn(4)]
    period = fourier_factor*scale
    coi=coi*dt*np.concatenate(([1.E-5],np.arange(1.,(n1+1.)/2.-1),np.flipud(np.arange(1,n1/2.)),[1.E-5])) # COI [Sec.3g]
    wave = wave[:,:n1]  # get rid of padding before returning
    return wave,period,scale,coi
    # end of code