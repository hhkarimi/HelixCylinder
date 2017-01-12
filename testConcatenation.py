# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 12:19:55 2016

Concatenate files for velocity data

@author: hkarimi
"""

#Test this with what c++ is reading in DER
import numpy as np

#%% Load data
#coords = np.loadtxt('./concatenated/coords.txt')
ux = np.loadtxt('./concatenated/ux.txt')
uy = np.loadtxt('./concatenated/uy.txt')
uz = np.loadtxt('./concatenated/uz.txt')
uxBar = ux; uyBar = uy; uzBar = uz;

#%%
Nr = 21; Nf = 41; Nz = 41; Nb = 21; sd = 3;
NperSD = Nr*Nf*Nz*Nb;
Ntot = Nr*Nf*Nz*Nb*sd;
#rArray = np.linspace(1.04e-4,1.0,11)
rArray = np.concatenate(([0.01], np.linspace(0.05, 1, Nr-1)), axis=0)
phiArray = np.linspace(0, 2*np.pi, Nf)
zArray = np.linspace(-1, 1, Nz)
bArray = np.linspace(0, 1, Nb)

def lowestindex (array, arrayLength, target):
    N = arrayLength
    index = -1;
    for i in range(0,N):
        if target < array[i]:
            index = i - 1
            break
    if index == -1:
        print("Error: target value is not within array.")
        print("target = {:g}, array[end] = {:g}".format(target,array[-1]))
    return index

rnorm = 0.525; phi = 4.1; z = 0.025; bnorm = 0.075; sd = 1;

# interpolate indices
nr = lowestindex( rArray, Nr, rnorm)
Dr = rArray[nr+1] - rArray[nr]
dr = rnorm - rArray[nr]

nf = lowestindex( phiArray, Nf, phi)
Df = phiArray[nf+1] - phiArray[nf]
df = phi - phiArray[nf]

nz = lowestindex( zArray, Nz, z)
Dz = zArray[nz+1] - zArray[nz]
dz = z - zArray[nz]

nb = lowestindex( bArray, Nb, bnorm)
Db = bArray[nb+1] - bArray[nb]
db = bnorm - bArray[nb]

# interpolation velocities
n = nr + nf*Nr + nz*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb
#np = (nr+1) + (nf+1)*Nr + (nz+1)*Nr*Nf + (nb+1)*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb
npr = (nr+1) + nf*Nr + nz*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
npf = nr + (nf+1)*Nr + nz*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
npz = nr + nf*Nr + (nz+1)*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
npb = nr + nf*Nr + nz*Nr*Nf + (nb+1)*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
uxInt = ux[n] + (ux[npr]-ux[n])/Dr * dr + (ux[npf]-ux[n])/Df * df + (ux[npz]-ux[n])/Dz * dz + (ux[npb]-ux[n])/Db * db;
uyInt = uy[n] + (uy[npr]-uy[n])/Dr * dr + (uy[npf]-uy[n])/Df * df + (uy[npz]-uy[n])/Dz * dz + (uy[npb]-uy[n])/Db * db;
uzInt = uz[n] + (uz[npr]-uz[n])/Dr * dr + (uz[npf]-uz[n])/Df * df + (uz[npz]-uz[n])/Dz * dz + (uz[npb]-uz[n])/Db * db;
print('ux,uy,uz[interpolationNaive] = {:g}, {:g}, {:g}'.format(uxInt,uyInt,uzInt))

# interpolate using multivariate linear interpolator
nr = lowestindex( rArray, Nr, rnorm)
Delta_r = rArray[nr+1] - rArray[nr];
dr = rnorm - rArray[nr];
Dr = rArray[nr+1] - rnorm;

nf = lowestindex( phiArray, Nf, phi);
Delta_f = phiArray[nf+1] - phiArray[nf];
df = phi - phiArray[nf];
Df = phiArray[nf+1] - phi;

nz = lowestindex( zArray, Nz, z);
Delta_z = zArray[nz+1] - zArray[nz];
dz = z - zArray[nz];
Dz = zArray[nz+1] - z;

nb = lowestindex( bArray, Nb, bnorm);
Delta_b = bArray[nb+1] - bArray[nb];
db = bnorm - bArray[nb];
Db = bArray[nb+1] - bnorm;
                
n = nr + nf*Nr + nz*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
npr = (nr+1) + nf*Nr + nz*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
npf = nr + (nf+1)*Nr + nz*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
npz = nr + nf*Nr + (nz+1)*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
npb = nr + nf*Nr + nz*Nr*Nf + (nb+1)*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
nprpf = (nr+1) + (nf+1)*Nr + nz*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
nprpz = (nr+1) + nf*Nr + (nz+1)*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
nprpb = (nr+1) + nf*Nr + nz*Nr*Nf + (nb+1)*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
npfpz = nr + (nf+1)*Nr + (nz+1)*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
npfpb = nr + (nf+1)*Nr + nz*Nr*Nf + (nb+1)*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
npzpb = nr + nf*Nr + (nz+1)*Nr*Nf + (nb+1)*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
nprpfpz = (nr+1) + (nf+1)*Nr + (nz+1)*Nr*Nf + nb*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
nprpfpb = (nr+1) + (nf+1)*Nr + nz*Nr*Nf + (nb+1)*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
nprpzpb = (nr+1) + nf*Nr + (nz+1)*Nr*Nf + (nb+1)*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
npfpzpb = nr + (nf+1)*Nr + (nz+1)*Nr*Nf + (nb+1)*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;
nprpfpzpb = (nr+1) + (nf+1)*Nr + (nz+1)*Nr*Nf + (nb+1)*Nr*Nf*Nz + sd*Nr*Nf*Nz*Nb;

deltaTotal = Delta_r * Delta_f * Delta_z * Delta_b;
uxInt = 1/deltaTotal * ( uxBar[n]*Dr*Df*Dz*Db  \
     + uxBar[npr] *dr*Df*Dz*Db + uxBar[npf] *Dr*df*Dz*Db  \
     + uxBar[npz] *Dr*Df*dz*Db + uxBar[npb] *Dr*Df*Dz*db  \
     + uxBar[nprpf] *dr*df*Dz*Db + uxBar[nprpz] *dr*Df*dz*Db  \
     + uxBar[nprpb] *dr*Df*Dz*db + uxBar[npfpz] *Dr*df*dz*Db  \
     + uxBar[npfpb] *Dr*df*Dz*db + uxBar[npzpb] *Dr*Df*dz*db  \
     + uxBar[nprpfpz] *dr*df*dz*Db + uxBar[nprpfpb] *dr*df*Dz*db  \
     + uxBar[nprpzpb] *dr*Df*dz*db + uxBar[npfpzpb] *Dr*df*dz*db  \
     + uxBar[nprpfpzpb] *dr*df*dz*db);
uyInt = 1/deltaTotal * ( uyBar[n]*Dr*Df*Dz*Db \
    + uyBar[npr] *dr*Df*Dz*Db + uyBar[npf] *Dr*df*Dz*Db  \
    + uyBar[npz] *Dr*Df*dz*Db + uyBar[npb] *Dr*Df*Dz*db  \
    + uyBar[nprpf] *dr*df*Dz*Db + uyBar[nprpz] *dr*Df*dz*Db  \
    + uyBar[nprpb] *dr*Df*Dz*db + uyBar[npfpz] *Dr*df*dz*Db  \
    + uyBar[npfpb] *Dr*df*Dz*db + uyBar[npzpb] *Dr*Df*dz*db  \
    + uyBar[nprpfpz] *dr*df*dz*Db + uyBar[nprpfpb] *dr*df*Dz*db  \
    + uyBar[nprpzpb] *dr*Df*dz*db + uyBar[npfpzpb] *Dr*df*dz*db  \
    + uyBar[nprpfpzpb] *dr*df*dz*db);
uzInt = 1/deltaTotal * ( uzBar[n]*Dr*Df*Dz*Db  \
    + uzBar[npr] *dr*Df*Dz*Db + uzBar[npf] *Dr*df*Dz*Db  \
    + uzBar[npz] *Dr*Df*dz*Db + uzBar[npb] *Dr*Df*Dz*db  \
    + uzBar[nprpf] *dr*df*Dz*Db + uzBar[nprpz] *dr*Df*dz*Db  \
    + uzBar[nprpb] *dr*Df*Dz*db + uzBar[npfpz] *Dr*df*dz*Db  \
    + uzBar[npfpb] *Dr*df*Dz*db + uzBar[npzpb] *Dr*Df*dz*db  \
    + uzBar[nprpfpz] *dr*df*dz*Db + uzBar[nprpfpb] *dr*df*Dz*db  \
    + uzBar[nprpzpb] *dr*Df*dz*db + uzBar[npfpzpb] *Dr*df*dz*db  \
    + uzBar[nprpfpzpb] *dr*df*dz*db);
    
print('ux,uy,uz[interpolationQuadLinear] = {:g}, {:g}, {:g}'.format(uxInt,uyInt,uzInt))
print('ux,uy,uz[{}] = {:g}, {:g}, {:g}'.format(n,ux[n],uy[n],uz[n]))
print('ux,uy,uz[{}] = {:g}, {:g}, {:g}'.format(nprpfpzpb,ux[nprpfpzpb],uy[nprpfpzpb],uz[nprpfpzpb]))

##### Use numpy interpolation packages
#coords_sp = (rArray, phiArray, zArray, bArray)
## convert ux_sd0 to a multidimensional array
#
#ux_sd0 = ux[ NperSD*(sd) : NperSD*(sd+1) ]
#coord = (rnorm,phi,z,bnorm)
#ux_sd0_sp = interpn(coords_sp,ux_sd0,coord)