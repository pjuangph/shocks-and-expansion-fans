import numpy as np 
import shockjump

def obshock(M:float,theta:float,gam:float=1.4,IsWeak:bool=True):
    """_summary_

    Args:
        M (float): _description_
        theta (float): _description_
        gam (float, optional): _description_. Defaults to 1.4.
        IsWeak (bool, optional): _description_. Defaults to True.
    """
    # input assignment and validation
	# Anderson, Modern Compressible Flow, pg 143
    if IsWeak:
        n=1
    else:
        n=0
    thetar = np.radians(theta)
    
    a = (M**2 - 1)
    b = 1 + M**2 * (gam - 1) / 2
    c = 1 + M**2 * (gam + 1) / 2
    lam = np.sqrt(a**2 - 3 * b * c * np.tan(thetar)**2) # Eq. 4.20
    x = (a**3 - 9*b*(1 + M**2*(gam-1)/2 + M**4*(gam+1)/4) * np.tan(thetar)**2) / lam**3
    tanB = (M**2 - 1 + 2*lam*np.cos( np.radians((720*n + np.arccos(x)*180/np.pi)/3))) / (3*b * np.tan(thetar))
    beta = np.degrees(np.arctan(tanB))
	
    assert not np.isnan(beta),'An oblique shock cannot exist under these conditions. Shock detached.'
	
    Mn1 = M * np.sin(np.radians(beta))
    
    Mn2 = shockjump.M2(Mn1,gam)
    M2 = Mn2 / np.sin(np.radians(beta - theta))
	
	## theta_max calculation
	# from Hady K. Joumaa, "Analytic Characterization of Oblique Shock Waves in Flows Around Wedges" (https://arxiv.org/pdf/1802.04763)
	
	# Analytical root of Eq. (5)
    theta_max = ((2**(1/2)*(2*M**8*gam - 24*M**4*gam - 32*M**2*gam + \
		(M**4*(gam + 1)*(8*M**2*gam + M**4*gam - 8*M**2 + M**4 + 16)**3)**(1/2) + \
        32*M**2 - 48*M**4 + 20*M**6 + M**8 - 8*M**4*gam**2 - 20*M**6*gam**2 + M**8*gam**2 - 32)**(1/2))/(2*(M**2*gam + M**2 + 2)**(3/2)*(M**2*gam - M**2 + 2)**(1/2)))
    theta_max = np.degrees(np.arctan(theta_max))
    
    if(theta == 0):
        M2 = 0
        Mn1 = 1
    return beta,Mn1,M2,theta_max

if __name__=='__main__':
    M=3.5 
    theta= 15
    gam=1.4
    beta,Mn1,M2,theta_max = obshock(M=3.5,theta=15,gam=1.4,IsWeak=True)
    print(f"Before Oblique Shock \n\t theta:{theta} M1:{M}")
    print(f"After Oblique Shock \n\t Beta:{beta} Mn1:{Mn1} M2:{M2} ThetaMax:{theta_max}")