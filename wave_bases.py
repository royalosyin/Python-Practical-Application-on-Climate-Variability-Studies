def wave_bases(mother,k,scale,param):
    """
    This is translation of wave_bases.m by Torrence and Gilbert P. Compo
 
    The folloing is the original README
 
%    WAVE_BASES  1D Wavelet functions Morlet, Paul, or DOG
%
%  [DAUGHTER,FOURIER_FACTOR,COI,DOFMIN] = ...
%      wave_bases(MOTHER,K,SCALE,PARAM);
%
%   Computes the wavelet function as a function of Fourier frequency,
%   used for the wavelet transform in Fourier space.
%   (This program is called automatically by WAVELET)
%
% INPUTS:
%
%    MOTHER = a string, equal to 'MORLET' or 'PAUL' or 'DOG'
%    K = a vector, the Fourier frequencies at which to calculate the wavelet
%    SCALE = a number, the wavelet scale
%    PARAM = the nondimensional parameter for the wavelet function
%
% OUTPUTS:
%
%    DAUGHTER = a vector, the wavelet function
%    FOURIER_FACTOR = the ratio of Fourier period to scale
%    COI = a number, the cone-of-influence size at the scale
%    DOFMIN = a number, degrees of freedom for each point in the wavelet power
%             (either 2 for Morlet and Paul, or 1 for the DOG)
%
%----------------------------------------------------------------------------
%   Copyright (C) 1995-1998, Christopher Torrence and Gilbert P. Compo
%   University of Colorado, Program in Atmospheric and Oceanic Sciences.
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.
%----------------------------------------------------------------------------
    """
    #import modules
    import numpy as np

    #
    mother = mother.upper()
    n = len(k)
    # define Heaviside step function
    def ksign(x):
        y=np.zeros_like(x)
        y[x>0]=1
        return y
    #
    if mother=='MORLET':  #-----------------------------------  Morlet
        if (param == -1): param = 6.
        k0 = param
        expnt = -(scale*k - k0)**2/2. *ksign(k)
        norm = np.sqrt(scale*k[1])*(np.pi**(-0.25))*np.sqrt(n)    # total energy=N   [Eqn(7)]
        daughter = norm*np.exp(expnt)
        daughter = daughter*ksign(k)  # Heaviside step function
        fourier_factor = (4.*np.pi)/(k0 + np.sqrt(2. + k0**2)) # Scale-->Fourier [Sec.3h]
        coi = fourier_factor/np.sqrt(2)            # Cone-of-influence [Sec.3g]
        dofmin = 2.                          # Degrees of freedom
    elif mother=='PAUL': #--------------------------------  Paul
        if (param == -1): param = 4.
        m = param
        expnt = -(scale*k)*ksign(k)
        norm = np.sqrt(scale*k[1])*(2.**m/np.sqrt(m*np.prod(np.arange(2,2*m))))*np.sqrt(n)
        daughter = norm*((scale*k)**m)*np.exp(expnt)
        daughter = daughter*ksign(k)      # Heaviside step function
        fourier_factor = 4*np.pi/(2.*m+1.)
        coi = fourier_factor*np.sqrt(2)
        dofmin = 2.
    elif mother=='DOG':  #--------------------------------  DOG
        if (param == -1): param = 2.
        m = param
        expnt = -(scale*k)**2 / 2.0
        from scipy.special import gamma 
        norm = np.sqrt(scale*k[1]/gamma(m+0.5))*np.sqrt(n)
        daughter = -norm*(1j**m)*((scale*k)**m)*np.exp(expnt);
        fourier_factor = 2.*np.pi*np.sqrt(2./(2.*m+1.))
        coi = fourier_factor/np.sqrt(2)
        dofmin = 1.
    else:
        raise Exception("Mother must be one of MORLET,PAUL,DOG")


    return daughter,fourier_factor,coi,dofmin 