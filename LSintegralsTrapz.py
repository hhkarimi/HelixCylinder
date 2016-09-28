#!/Users/hkarimi/anaconda/lib/python3.5

##################################### LOAD MODULES #####################################
    # for performance timing
import time # use as: tic = time.time(); elapsedTime = time.time() - tic
    # numerical modules
import numpy as np
from numpy import cos,sin,sqrt
from scipy import integrate,special
#####################################
import matplotlib.pyplot as plt # main plotting module
plt.rc('text',usetex=True) # for latex
plt.rc('font',family='serif') # for latex
    # for sending data to MATLAB
import scipy.io as sio

# treating division by 0
np.seterr(divide='ignore')

################  use analytical solution to stokeslet flow in a pipe ################

################  function definitions ################
def f(veldir,stokesdir,l,R,k):
    lR = l*R
    lb = l*b
    # modified bessel of the first kind
    Ik_lb = special.iv(k,lb)
    Ikm1_lb = special.iv(k-1,lb)
    Ikp1_lb = special.iv(k+1,lb)
    Ikm2_lb = special.iv(k-2,lb)
    Ikp2_lb = special.iv(k+2,lb)
        # derivatives
    dIk_lb = special.ivp(k,lb,1)
    dIkm1_lb = special.ivp(k-1,lb,1)
    dIkp1_lb = special.ivp(k+1,lb,1)
    # modified bessel of the second kind
    Kk_lb = special.kn(k,lb)
    Kk_lR = special.kn(k,lR)
    Kkm1_lR = special.kn(k-1,lR)
    Kkp1_lR = special.kn(k+1,lR)
        # derivatives
    dKk_lR = special.kvp(k,lR,1)

    if R >= b:
        if (stokesdir==1) and (veldir==1):
            return Kk_lR*(Ik_lb-lb*dIk_lb) + lR/2*(Kkm1_lR*dIkm1_lb+Kkp1_lR*dIkp1_lb)
        elif ( stokesdir==1 and veldir==2 ) or ( stokesdir==2 and veldir==1 ):
            return 1/4*lR*(Kkm1_lR*Ikm2_lb-Kkp1_lR*Ikp2_lb) - 1/2*k*Kk_lR*Ik_lb
        elif ( stokesdir==1 and veldir==3 ) or ( stokesdir==3 and veldir==1 ):
            return 1/2*lR*(Kkm1_lR*Ikm1_lb+Kkp1_lR*Ikp1_lb) - lb*Kk_lR*Ik_lb
        elif (stokesdir==2) and (veldir==2):
            return 1/2*Ik_lb*(2*Kk_lR-lR*dKk_lR) - 1/4*lR*(Kkp1_lR*Ikp2_lb+Kkm1_lR*Ikm2_lb)
        elif ( stokesdir==2 and veldir==3 ) or ( stokesdir==3 and veldir==2 ):
            return 1/2*lR*(Kkm1_lR*Ikm1_lb-Kkp1_lR*Ikp1_lb)
        elif (stokesdir==3) and (veldir==3):
            return Ik_lb*(2*Kk_lR+lR*dKk_lR) + lb*dIk_lb*Kk_lR
    else: # if R < b, switch K and I or witch b and R?
        # modified bessel of the first kind
        Ik_lb = special.iv(k,lb)
        Ik_lR = special.iv(k,lR)
        Ikm1_lb = special.iv(k-1,lb)
        Ikm1_lR = special.iv(k-1,lR)
        Ikp1_lb = special.iv(k+1,lb)
        Ikp1_lR = special.iv(k+1,lR)
        Ikm2_lb = special.iv(k-2,lb)
        Ikm2_lR = special.iv(k-2,lR)
        Ikp2_lb = special.iv(k+2,lb)
        Ikp2_lR = special.iv(k+2,lR)
            # derivatives
        dIk_lb = special.ivp(k,lb,1)
        dIk_lR = special.ivp(k,lR,1)
        dIkm1_lb = special.ivp(k-1,lb,1)
        dIkp1_lb = special.ivp(k+1,lb,1)
        # modified bessel of the second kind
        Kk_lb = special.kn(k,lb)
        Kk_lR = special.kn(k,lR)
        Kkm1_lR = special.kn(k-1,lR)
        Kkm1_lb = special.kn(k-1,lb)
        Kkp1_lR = special.kn(k+1,lR)
        Kkp1_lb = special.kn(k+1,lb)
        Kkm2_lb = special.kn(k-2,lb)
        Kkp2_lb = special.kn(k+2,lb)
            # derivatives
        dKk_lR = special.kvp(k,lR,1)
        dKk_lb = special.kvp(k,lb,1)
        dKkm1_lb = special.kvp(k-1,lb,1)
        dKkp1_lb = special.kvp(k+1,lb,1)
        if (stokesdir==1) and (veldir==1):
            return Ik_lR*(Kk_lb-lb*dKk_lb) + lR/2*(Ikm1_lR*dKkm1_lb+Ikp1_lR*dKkp1_lb)
        elif ( stokesdir==1 and veldir==2 ) or ( stokesdir==2 and veldir==1 ):
            return -1/4*lR*(Ikm1_lR*Kkm2_lb-Ikp1_lR*Kkp2_lb) - 1/2*k*Ik_lR*Kk_lb
        elif ( stokesdir==1 and veldir==3 ) or ( stokesdir==3 and veldir==1 ):
            return ( 1/2*lR*(Ikm1_lR*Kkm1_lb+Ikp1_lR*Kkp1_lb) - lb*Ik_lR*Kk_lb )
        elif (stokesdir==2) and (veldir==2):
            return 1/2*Kk_lb*(2*Ik_lR-lR*dIk_lR) + 1/4*lR*(Ikp1_lR*Kkp2_lb+Ikm1_lR*Kkm2_lb)
        elif ( stokesdir==2 and veldir==3 ) or ( stokesdir==3 and veldir==2 ):
            return 1/2*lR*(Ikm1_lR*Kkm1_lb-Ikp1_lR*Kkp1_lb)
        elif (stokesdir==3) and (veldir==3):
            return Kk_lb*(2*Ik_lR+lR*dIk_lR) + lb*dKk_lb*Ik_lR

def psi(stokesdir,k,l):
    s = l*R0
    
    # modified bessel of the first kind
    Ik_s = special.iv(k,s)
    Ikm1_s = special.iv(k-1,s)
    Ikp1_s = special.iv(k+1,s)
        # derivatives
    dIk_s = special.ivp(k,s,1)
    dIkm1_s = special.ivp(k-1,s,1)
    dIkp1_s = special.ivp(k+1,s,1)
    
    Delta = l*(s*Ik_s*( dIkm1_s*Ikp1_s+Ikm1_s*dIkp1_s) - 2*(Ik_s+s*dIk_s)*Ikm1_s*Ikp1_s )
    return ( -s*(dIkm1_s*Ikp1_s+Ikm1_s*dIkp1_s)*L(stokesdir,k,l)
        - (Ik_s+s*dIk_s)*(Ikm1_s*H(stokesdir,k,l)+Ikp1_s*G(stokesdir,k,l)) ) / Delta

def omega(stokesdir,k,l):
    s_omega = -1 if stokesdir == 2 else 1
    s = l*R0

    # modified bessel of the first kind
    Ik_s = special.iv(k,s)
    Ikm1_s = special.iv(k-1,s)
    Ikp1_s = special.iv(k+1,s)
        # derivatives
    dIk_s = special.ivp(k,s,1)
    dIkm1_s = special.ivp(k-1,s,1)
    dIkp1_s = special.ivp(k+1,s,1)
    
    Delta = l*(s*Ik_s*( dIkm1_s*Ikp1_s+Ikm1_s*dIkp1_s) - 2*(Ik_s+s*dIk_s)*Ikm1_s*Ikp1_s )
    return s_omega*( (s*Ik_s*dIkp1_s-Ikp1_s*(Ik_s+s*dIk_s))*G(stokesdir,k,l) + ((Ik_s+s*dIk_s)*Ikm1_s-s*Ik_s*dIkm1_s)*H(stokesdir,k,l) + (Ikm1_s*dIkp1_s-dIkm1_s*Ikp1_s)*s*L(stokesdir,k,l) ) / Delta

def pi(stokesdir,k,l):
    s = l*R0
    
    # modified bessel of the first kind
    Ik_s = special.iv(k,s)
    Ikm1_s = special.iv(k-1,s)
    Ikp1_s = special.iv(k+1,s)
        # derivatives
    dIk_s = special.ivp(k,s,1)
    dIkm1_s = special.ivp(k-1,s,1)
    dIkp1_s = special.ivp(k+1,s,1)
    
    Delta = l*(s*Ik_s*( dIkm1_s*Ikp1_s+Ikm1_s*dIkp1_s) - 2*(Ik_s+s*dIk_s)*Ikm1_s*Ikp1_s )
    return ( Ik_s*(Ikp1_s*G(stokesdir,k,l)+Ikm1_s*H(stokesdir,k,l)) + 2*L(stokesdir,k,l)*Ikm1_s*Ikp1_s ) / Delta

def H(stokesdir,k,l):
    if stokesdir == 1: # x
        return -1/(4*np.pi**2*mu) * ( f(1,1,l,R0,k+1)+f(2,1,l,R0,k+1) )
    elif stokesdir == 2: # y
        return -1/(4*np.pi**2*mu) * ( f(1,2,l,R0,k+1)-f(2,2,l,R0,k+1) )
    elif stokesdir == 3: # z
        return -1/(4*np.pi**2*mu) * ( f(1,3,l,R0,k+1)+f(2,3,l,R0,k+1) )

def G(stokesdir,k,l):
    if stokesdir == 1: # x
        return -1/(4*np.pi**2*mu) * ( f(1,1,l,R0,k-1)-f(2,1,l,R0,k-1) )
    elif stokesdir == 2: # y
        return -1/(4*np.pi**2*mu) * ( f(1,2,l,R0,k-1)+f(2,2,l,R0,k-1) )
    elif stokesdir == 3: # z
        return -1/(4*np.pi**2*mu) * ( f(1,3,l,R0,k-1)-f(2,3,l,R0,k-1) )

def L(stokesdir,k,l):
    if stokesdir == 1:
        return -1/(4*np.pi**2*mu) * f(3,1,l,R0,k)
    elif stokesdir == 2:
        return -1/(4*np.pi**2*mu) * f(3,2,l,R0,k)
    elif stokesdir == 3:
        return 1/(4*np.pi**2*mu) * f(3,3,l,R0,k)

def gR_int(stokesdir,k,Phi,l,R):
    j = stokesdir
    Ik = special.iv(k,l*R) # modified bessel of the first kind and derivatives
    dIk = special.ivp(k,l*R,1)
    ddIk = special.ivp(k,l*R,2)
    if stokesdir == 1:
        tmp = cos(k*Phi)*psi(j,k,l)*l*dIk*cos(l*z) - sin(k*Phi-np.pi/2)*omega(j,k,l)*(k/R)*Ik*cos(l*z) + cos(k*Phi)*(l**2)*R*pi(j,k,l)*ddIk*cos(l*z)
        return tmp
    elif stokesdir == 2:
        tmp = cos(k*Phi-np.pi/2)*psi(j,k,l)*l*dIk*cos(l*z) - sin(k*Phi)*omega(j,k,l)*(k/R)*Ik*cos(l*z) + cos(k*Phi-np.pi/2)*l**2*R*pi(j,k,l)*ddIk*cos(l*z)
        return tmp
    elif stokesdir == 3:
        tmp = cos(k*Phi)*psi(j,k,l)*l*dIk*cos(l*z-np.pi/2) - sin(k*Phi-np.pi/2)*omega(j,k,l)*(k/R)*Ik*cos(l*z-np.pi/2) + cos(k*Phi)*(l**2)*R*pi(j,k,l)*ddIk*cos(l*z-np.pi/2)
        return tmp

def gPhi_int(stokesdir,k,Phi,l,R):
    j = stokesdir
    Ik = special.iv(k,l*R) # modified bessel of the first kind and derivatives
    dIk = special.ivp(k,l*R,1)
    ddIk = special.ivp(k,l*R,2)
    if stokesdir == 1:
        tmp = (-sin(k*Phi)*k/R*psi(j,k,l)*Ik*cos(l*z)  - cos(k*Phi-np.pi/2)*omega(j,k,l)*l*dIk*cos(l*z)                     + sin(k*Phi)*k/R*pi(j,k,l)*Ik*cos(l*z) -sin(k*Phi)*k*pi(j,k,l)*l*dIk*cos(l*z) )
        return tmp
    elif stokesdir == 2:
        tmp = -sin(k*Phi-np.pi/2)*k/R*psi(j,k,l)*Ik*cos(l*z) - cos(k*Phi)*omega(j,k,l)*l*dIk*cos(l*z) + sin(k*Phi-np.pi/2)*k/R*pi(j,k,l)*Ik*cos(l*z) - sin(k*Phi-np.pi/2)*k*pi(j,k,l)*l*dIk*cos(l*z)
        return tmp
    elif stokesdir == 3:
        tmp = -sin(k*Phi)*k/R*psi(j,k,l)*Ik*cos(l*z-np.pi/2) - cos(k*Phi-np.pi/2)*omega(j,k,l)*l*dIk*cos(l*z-np.pi/2)  + sin(k*Phi)*k/R*pi(j,k,l)*Ik*cos(l*z-np.pi/2) - sin(k*Phi)*k*pi(j,k,l)*l*dIk*cos(l*z-np.pi/2)
        return tmp

def gz_int(stokesdir,k,Phi,l,R):
    j = stokesdir
    Ik = special.iv(k,l*R) # modified bessel of the first kind and derivatives
    dIk = special.ivp(k,l*R,1)
    ddIk = special.ivp(k,l*R,2)
    if stokesdir == 1:
        tmp = -cos(k*Phi)*l*psi(j,k,l)*Ik*sin(l*z) - cos(k*Phi)*(l**2*R*dIk+l*Ik)*pi(j,k,l)*sin(l*z)
        return tmp
    elif stokesdir == 2:
        tmp = -cos(k*Phi-np.pi/2)*l*psi(j,k,l)*Ik*sin(l*z) - cos(k*Phi-np.pi/2)*(l**2*R*dIk+l*Ik)*pi(j,k,l)*sin(l*z)
        return tmp
    elif stokesdir == 3:
        tmp = -cos(k*Phi)*l*psi(j,k,l)*Ik*sin(l*z-np.pi/2) - cos(k*Phi)*(l**2*R*dIk+l*Ik)*pi(j,k,l)*sin(l*z-np.pi/2)
        return tmp

##### Stokeslet position and evaluation point
R0 = 1 # radius of cylinder
mu = 1 # viscosity of fluid
b = 0.5 # eccentricity of stokeslet
##### evaluation points
Reval = np.linspace(0.001,R0,15) # radial position from stokeslet
Phi = 37*np.pi/180 # angular displacement from stokeslet
z = 0.4 # axial position from stokeslet
fx = 0.0; fy = 1.0; fz = 0.0;

PhiDeg = Phi * 180/np.pi
print("(z,Phi = {:g},{:g}".format(z,PhiDeg) + u'\u00B0)')

#### loop over all force directions
tic = time.time();
uTotNorm = np.array([])
uTot = np.array([]).reshape(3,0)
vTot = np.array([]).reshape(3,0)
gTot = np.array([]).reshape(3,0)
K = 5 # summation index K -> infinity
l = np.linspace(0.0001, 30/R0, num=250)
for nr in range(0,len(Reval)):
    R = Reval[nr]
    # initialize velocities for summations
    vx = 0; vy = 0; vz = 0;
    gR = 0; gPhi = 0; gz = 0;
    for j in range(1,4): # loop over all 3 stokeslet directions
        for k in range(-K,K+1):
            if j == 1 and fx != 0: # due to stokeslet in x direction
                vx = vx + fx*cos(k*Phi) * integrate.trapz( f(1,j,l,R,k)*cos(l*z) , l)
                vy = vy + fx*sin(k*Phi) * integrate.trapz( f(2,j,l,R,k)*cos(l*z) , l)
                vz = vz + fx*cos(k*Phi) * integrate.trapz( f(3,j,l,R,k)*sin(l*z) , l)
                gR = gR + fx*integrate.trapz( gR_int(j,k,Phi,l,R), l)
                gPhi = gPhi + fx*integrate.trapz( gPhi_int(j,k,Phi,l,R), l)
                gz = gz + fx*integrate.trapz( gz_int(j,k,Phi,l,R), l)
                if any( np.isnan([vx,vy,vz])): sys.exit("error in v")
                if any( np.isnan([gR,gPhi,gz])):
                    print("gR,gPhi,gz = {:g},{:g},{:g}".format(gR,gPhi,gz))
                    sys.exit("error in g")
            elif j == 2 and fy != 0: # due to stokeslet in y direction
                vx = vx + fy*sin(k*Phi) * integrate.trapz( f(1,j,l,R,k)*cos(l*z) , l)
                vy = vy + fy*cos(k*Phi) * integrate.trapz( f(2,j,l,R,k)*cos(l*z) , l)
                vz = vz + fy*sin(k*Phi) * integrate.trapz( f(3,j,l,R,k)*sin(l*z) , l)
                gR = gR + fy*integrate.trapz( gR_int(j,k,Phi,l,R), l)
                gPhi = gPhi + fy*integrate.trapz( gPhi_int(j,k,Phi,l,R), l)
                gz = gz + fy*integrate.trapz( gz_int(j,k,Phi,l,R), l)
            elif j == 3 and fz !=0: # due to stokeslet in z direction
                vx = vx + fz*cos(k*Phi) * integrate.trapz( f(1,j,l,R,k)*sin(l*z) , l)
                vy = vy + fz*sin(k*Phi) * integrate.trapz( f(2,j,l,R,k)*sin(l*z) , l)
                vz = vz + fz*cos(k*Phi) * integrate.trapz( f(3,j,l,R,k)*cos(l*z) , l)
                gR = gR + fz*integrate.trapz( gR_int(j,k,Phi,l,R), l)
                gPhi = gPhi + fz*integrate.trapz( gPhi_int(j,k,Phi,l,R), l)
                gz = gz + fz*integrate.trapz( gz_int(j,k,Phi,l,R), l)
        # END k summation
        
        gx = gR*cos(Phi)-gPhi*sin(Phi)
        gy = gR*sin(Phi)+gPhi*cos(Phi)
        gz = gz
    # END loop over j (stokesdir) direction
    
    v = np.array([ [vx],[vy],[vz] ])/(4*np.pi**2*mu)
    g = np.array([ [gx],[gy],[gz] ])
    u = v + g
    uTot = np.append( uTot, u , axis=1 )
    vTot = np.append( vTot, v , axis=1)
    gTot = np.append( gTot, g , axis=1 )
    uTotNorm = np.append( uTotNorm, np.linalg.norm(u) )
    print(uTotNorm[nr])
elapsedTime = time.time() - tic
print("elapsed time of computation = {:3f}".format(elapsedTime))


########################### plotting the solution ##########################
plt.plot(uTotNorm,Reval)
plt.xlabel(r'$|u|$')
plt.ylabel(r'$R$')
#plt.show()

fig2 = plt.figure()
xLabels = ['$u_x$','$u_y$','$u_z$']
for n in range(0,3):
    ax = fig2.add_subplot(1,3,n+1)
    ax.plot(8*mu*uTot[n][:],Reval,'k')
#    ax.plot(vTot[n][:],Reval,'r')
#    ax.plot(gTot[n][:],Reval,'b')
#    if n == 1: ax.legend([r'$u_\textrm{tot}$',r'$v_\textrm{tot}$',r'$g_\textrm{tot}$'])
    ax.plot( [0,0] , [0,R0] ,'--g' )
    ax.set_xlabel(xLabels[n])
    if n == 0: ax.set_ylabel(r'$R$')
    ax.locator_params(axis='x',nbins=5)
fig2.set_tight_layout

plt.show()

figureName = eval("'phi' + '{:d}'.format(int(PhiDeg)) + '.pdf'")
fig2.savefig(figureName,format='pdf')

#### SAVE DATA AS FLOAT32!!! 
# u[0][:].astype('float32')
