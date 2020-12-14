import numpy as np

def streamfunc(a,b,Pa,Pb,thetatp,z=0):
    f = 1e-4
    theta00 = 300
    g = 10
    Ns = 2e-2
    Nt = 1e-2
    pâté = g*(Ns-Nt)/(theta00*Ns*Nt)
    
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

    psi = np.fft.ifft2(psihat).real
    
    return psi

def geostwind(a,b,Pa,Pb,psi):
    ug = -(np.roll(psi,-1,0)-np.roll(psi,1,0))/(2*a/Pa) # centered finite differences
    vg = (np.roll(psi,-1,1)-np.roll(psi,1,1))/(2*b/Pb)
    
    return ug,vg