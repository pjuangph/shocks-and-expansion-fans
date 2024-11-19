from __future__ import absolute_import
from typing import Union
import numpy.typing as npt 
import numpy as np 
from scipy.optimize import minimize_scalar

def P2_P1(M1:Union[npt.NDArray,float], M2:Union[npt.NDArray,float], gamma:float=1.4) -> Union[npt.NDArray,float]:
    """Returns P2/P1 ratio for Rayleigh Flow

    Args:
        M1 (Union[npt.NDArray,float]): Mach 1 
        M2 (Union[npt.NDArray,float]): Mach 2 
        gamma (float, optional): ratio of specific heats. Defaults to 1.4.

    Returns:
        Union[npt.NDArray,float]: Returns P2_P2 Ratio 
    """
    return ((1 + gamma*M1**2) / (1 + gamma*M2**2))

def T2_T1(M1:Union[npt.NDArray,float], M2:Union[npt.NDArray,float], gamma:float=1.4) ->Union[npt.NDArray,float]:
    """Return T2/T1 ratio for rayleigh flow

    Args:
        M1 (Union[npt.NDArray,float]): Incoming mach number
        M2 (Union[npt.NDArray,float]): Exit mach number 
        gamma (float, optional): Ratio of specific heats. Defaults to 1.4.

    Returns:
        Union[npt.NDArray,float]: Returns T2/T1 
    """
    return (
        ((1 + gamma*M1**2) / (1 + gamma*M2**2))**2 *
        (M2**2 / M1**2)
        )
            
def rho2_rho1(M1:Union[npt.NDArray,float], M2:Union[npt.NDArray,float], gamma:float=1.4) -> Union[npt.NDArray,float]:
    """Return rho2 rho1 ratio for rayleigh flow

    Args:
        M1 (Union[npt.NDArray,float]): Incoming mach number
        M2 (Union[npt.NDArray,float]): Exit mach number 
        gamma (float, optional): Ratio of specific heats. Defaults to 1.4.

    Returns:
        Union[npt.NDArray,float]: density ratio rho2 rho1
    """
    return (
        ((1 + gamma*M2**2) / (1 + gamma*M1**2)) *
        (M1**2 / M2**2)
        )
            
def T02_T01(M1:Union[npt.NDArray,float], M2:Union[npt.NDArray,float], gamma:float=1.4) -> Union[npt.NDArray,float]:
    '''Return Tt2/Tt1 for Rayleigh flow'''
    return (
        ((1 + gamma*M1**2) / (1 + gamma*M2**2))**2 * 
            (M2 / M1)**2 * (
            (1 + 0.5*(gamma-1)*M2**2) /
            (1 + 0.5*(gamma-1)*M1**2)
            )
        )

def P02_P01(M1:Union[npt.NDArray,float], M2:Union[npt.NDArray,float], gamma:float=1.4) -> Union[npt.NDArray,float]:
    '''Return pt2/pt1 for Rayleigh flow'''
    return (
        ((1 + gamma*M1**2) / (1 + gamma*M2**2)) * (
            (1 + 0.5*(gamma-1)*M2**2) /
            (1 + 0.5*(gamma-1)*M1**2)
            )**(gamma / (gamma - 1))
        )

def P_P_sonic(mach:float, gamma:float=1.4):
    '''Return p/p* for Rayleigh flow'''
    return ((1 + gamma) / (1 + gamma*mach**2))

def T_T_sonic(mach, gamma:float=1.4):
    '''Return T/T* for Rayleigh flow'''
    return (
        mach**2 * (1 + gamma)**2 / 
        (1 + gamma*mach**2)**2
        )
            
def rho_rho_sonic(mach:float, gamma:float=1.4):
    '''Return rho/rho* for Rayleigh flow'''
    return ((1 + gamma*mach**2) / ((1 + gamma) * mach**2))
            
def T0_T0_sonic(mach:float, gamma:float=1.4):
    '''Return Tt/Tt* for Rayleigh flow'''
    return (
        2*(1 + gamma)*mach**2 * 
        (1 + 0.5*(gamma - 1)*mach**2) / (1 + gamma*mach**2)**2
        )

def P0_P0_sonic(mach:float, gamma:float=1.4):
    '''Return pt/pt* for Rayleigh flow'''
    return (
        (1 + gamma) * (
            (1 + 0.5*(gamma-1)*mach**2) / (0.5*(gamma+1))
            )**(gamma / (gamma - 1)) /
            (1 + gamma*mach**2)
        )

def Mach_T02_T01(T02T01:float,M1:float,gamma:float=1.4,IsSupersonic:bool=True):
    """Find the mach number for a given T02/T01 ratio and mach number 

    Args:
        T02T01 (float): Ratio of T02/T01
        M1 (float): Incoming mach number 
        gamma (float, optional): Ratio of specific heats. Defaults to 1.4.
        IsSupersonic (bool, optional): Use supersonic solution. Defaults to True.
    """
    def f(M2):
        return np.abs(T02_T01(M1, M2, gamma) - T02T01)
    if IsSupersonic:
        res = minimize_scalar(f,bounds=[1.01,5])
    else:
        res = minimize_scalar(f,bounds=[0.01,0.99])
    return res.x 
    