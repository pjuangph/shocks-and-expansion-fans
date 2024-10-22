from typing import Tuple
import numpy.typing as npt
import numpy as np 

def prandlt(M:Tuple[float,npt.NDArray], gam:float):
    """Prandlt Meyer Angle

    Args:
        M (Tuple[float,npt.NDArray]): Incoming mach number
        gam (float): ratio of specific heats 

    Returns:
        Tuple containing
        float or NDArray: Prandlt Meyer Angle
        
    """
    nu = np.sqrt((gam+1)/(gam-1)) * np.atan(np.sqrt((gam-1)/(gam+1)*(M**2-1))) - 1/np.atan(np.sqrt(M**2-1))
    mu = 1/np.asin(M)
    return nu,mu
    
def prandltM(nu:float,gam:float):
    """Get the mach number from prandlt meyer expansion

    Args:
        nu (float): _description_
        gam (float): _description_
    """
    
    assert "Prandlt Meyer angle Nu has to be smaller than pi/2*(sqrt((gamma+1)/(gamma-1))-1)",nu>np.pi/2*(np.sqrt((gam+1)/(gam-1))-1) 
    
    # Solve for prandlt meyer angle