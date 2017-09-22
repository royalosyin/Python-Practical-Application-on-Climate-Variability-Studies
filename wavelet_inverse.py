def wavelet_inverse(wave, scale, dt, dj=0.25, mother="MORLET",param=-1):
    """Inverse continuous wavelet transform
    Torrence and Compo (1998), eq. (11)

    INPUTS
        waves (array like):
          WAVE is the WAVELET transform. This is a complex array.
          
        scale (array like):
           the vector of scale indices 
        dt (float) :
            amount of time between each original value, i.e. the sampling time.
        dj (float, optional) :
            the spacing between discrete scales. Default is 0.25.
           A smaller # will give better scale resolution, but be slower to plot.
        mother (string, optional) :
            the mother wavelet function.
             The choices are 'MORLET', 'PAUL', or 'DOG'
         PARAM = the mother wavelet parameter.
            For 'MORLET' this is k0 (wavenumber), default is 6.
            For 'PAUL' this is m (order), default is 4.
            For 'DOG' this is m (m-th derivative), default is 2.    

    OUTPUTS
        iwave (array like) :
            Inverse wavelet transform.
    """
    import numpy as np
    
    j1, n = wave.shape
    J1 = len(scale)
    if not j1 == J1:
        print j1,n,J1
        raise Exception("Input array are inconsistent")
    sj = np.dot(scale.reshape(len(scale),1),np.ones((1,n)))
    #
    mother = mother.upper()
    
    # psi0 comes from Table 1,2 Torrence and Compo (1998)
    # Cdelta comes from Table 2 Torrence and Compo (1998)
    if mother=='MORLET':  #-----------------------------------  Morlet
        if (param == -1): param = 6.
        psi0=np.pi**(-0.25)
        if param==6.:
            Cdelta = 0.776
    elif mother=='PAUL': #--------------------------------  Paul
        if (param == -1): param = 4.
        m = param   
        psi0=np.real(2.**m*1j**m*np.prod(np.arange(2, m + 1))/np.sqrt(np.pi*np.prod(np.arange(2,2*m+1)))*(1**(-(m+1))))
        if m==4.:
           Cdelta = 1.132 
    elif mother=='DOG':  #--------------------------------  DOG
        if (param == -1): param = 2.
        m = param
        from scipy.special import gamma 
        from numpy.lib.polynomial import polyval
        from scipy.special.orthogonal import hermitenorm
        p = hermitenorm(m)
        psi0=(-1)**(m+1)/np.sqrt(gamma(m+0.5))*polyval(p, 0)
        print psi0
        if m==2.:
            Cdelta=3.541
        if m==6.:
            Cdelta=1.966
    else:
        raise Exception("Mother must be one of MORLET,PAUL,DOG")
    
    #eq. (11) in Torrence and Compo (1998)
    iwave = dj * np.sqrt(dt) / Cdelta /psi0 * (np.real(wave) / sj**0.5).sum(axis=0) 
    return iwave