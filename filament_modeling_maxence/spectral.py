import numpy as np

def geostwind(a,b,thetatp,z=0,fourier=False):
    f = 1e-4
    theta00 = 300
    g = 10
    Ns = 2e-2
    Nt = 1e-2
    pâté = g*(Ns-Nt)/(theta00*Ns*Nt)
    
    Pa, Pb = thetatp.shape
    
    thetatphat = np.fft.fft2(thetatp)
        
    freqx = np.fft.fftfreq(Pa, a/Pa)
    freqy = np.fft.fftfreq(Pb, b/Pb)
    wavevect = np.array([[(2*np.pi*m,2*np.pi*n) for n in freqy] for m in freqx])
    
    N = Ns if z>0 else Nt
    psihat = thetatphat
    for m in range(Pa):
        for n in range(Pb):
            K = np.linalg.norm(wavevect[m,n])
            if K==0:
                psihat[m,n] = 0
            else:
                psihat[m,n] *= pâté/K*np.exp(-N*K/f*np.abs(z))

    if fourier:
        ug = -np.fft.ifft2(psihat*1j*wavevect[:,:,1]).real
        vg = np.fft.ifft2(psihat*1j*wavevect[:,:,0]).real
    
    else:
        psi = np.fft.ifft2(psihat).real
        ug = -(np.roll(psi,-1,1)-np.roll(psi,1,1))/(2*a/Pa)
        vg = (np.roll(psi,-1,0)-np.roll(psi,1,0))/(2*b/Pb)
    
    return ug,vg

def vertwind(a,b,thetatp,thetatpprev,dt,z=0):
    f = 1e-4
    theta00 = 300
    g = 10
    Ns = 2e-2
    Nt = 1e-2
    pâté = g*(Ns-Nt)/(theta00*Ns*Nt)
    
    Pa, Pb = thetatp.shape
    
    thetatphat = np.fft.fft2(thetatp)
    thetatpprevhat = np.fft.fft2(thetatpprev)
        
    freqx = np.fft.fftfreq(Pa, a/Pa)
    freqy = np.fft.fftfreq(Pb, b/Pb)
    wavevect = np.array([[(2*np.pi*m,2*np.pi*n) for n in freqy] for m in freqx])
    
    N = Ns if z>0 else Nt
    psihat = thetatphat
    psiprevhat = thetatpprevhat
    thetazhat = np.zeros(thetatphat.shape,dtype=complex)
    thetazprevhat = np.zeros(thetatpprevhat.shape,dtype=complex)
    for m in range(Pa):
        for n in range(Pb):
            K = np.linalg.norm(wavevect[m,n])
            if K==0:
                psihat[m,n] = 0
                psiprevhat[m,n] = 0
                thetazhat[m,n] = 0
                thetazprevhat[m,n] = 0
            else:
                psihat[m,n] *= pâté/K*np.exp(-N*K/f*np.abs(z))
                psiprevhat[m,n] *= pâté/K*np.exp(-N*K/f*np.abs(z))
                thetazhat[m,n] = theta00/g*(-np.sign(z))*N*K*psihat[m,n]
                thetazprevhat[m,n] = theta00/g*(-np.sign(z))*N*K*psiprevhat[m,n]

    ug = -np.fft.ifft2(psihat*1j*wavevect[:,:,1]).real
    vg = np.fft.ifft2(psihat*1j*wavevect[:,:,0]).real
    
    thetaz = np.fft.ifft2(thetazhat).real
    thetazprev = np.fft.ifft2(thetazprevhat).real
    dtthetaz = (thetaz-thetazprev)/dt
    dxthetaz = np.fft.ifft2(thetazhat*1j*wavevect[:,:,0]).real
    dythetaz = np.fft.ifft2(thetazhat*1j*wavevect[:,:,1]).real
    
    w = -dtthetaz - ug*dxthetaz - vg*dythetaz
    w *= g/(N**2 * theta00)
    
    return w