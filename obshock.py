from scipy.optimize import minimize_scalar
import numpy as np 

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
        
    a = (M**2 - 1)
    b = 1 + M**2 * (gam - 1) / 2
    c = 1 + M**2 * (gam + 1) / 2
    lam = np.sqrt(a**2 - 3 * b * c * np.tan(theta)**2) #  Eq. 4.20
    x = (a**3 - 9*b*(1 + M**2*(gam-1)/2 + M**4*(gam+1)/4) * np.tan(theta)**2) / lam**3
    tanB = (M**2 - 1 + 2*lam*np.cos((4*np.pi*n + np.acos(x)) / 3)) / (3*b * np.atan(theta))
    beta = np.atan(tanB)
	
    assert np.imag(beta) ==0,'An oblique shock cannot exist under these conditions. Shock detached.'
	
    Mn1 = M * np.sin(beta)
    ~,~,~,~,Mn2 = flownormalshock(gam, Mn1)
    M2 = Mn2 / np.sin(beta - theta)
	
	## theta_max calculation
	# from Hady K. Joumaa, "Analytic Characterization of Oblique Shock Waves in Flows Around Wedges" (https://arxiv.org/pdf/1802.04763)
	
	# Analytical root of Eq. (5)
	theta_max = atand((2^(1/2)*(2*M.^8*gam - 24*M.^4*gam - 32*M.^2*gam +...
		(M.^4*(gam + 1).*(8*M.^2*gam + M.^4*gam - 8*M.^2 + M.^4 + 16).^3).^(1/2) +...
		32*M.^2 - 48*M.^4 + 20*M.^6 + M.^8 - 8*M.^4*gam^2 - 20*M.^6*gam^2 +...
		M.^8*gam^2 - 32).^(1/2))./(2*(M.^2*gam + M.^2 + 2).^(3/2).*(M.^2*gam - M.^2 + 2).^(1/2)));
	if(theta == 0)
		M2 = 0;
		Mn1 = 1;
	end