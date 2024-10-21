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
    pass