import numpy as np
from numpy import log, exp #Frequently used math function
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import minimize,shgo,curve_fit,Bounds

#Constants & Math Functions
π=np.pi

def log(a):
    if a>0:
        value=np.log(a)
    else:
        value=0.  # In the project, log function always appears in the form a log(a)
    return value
#Operators

##1st derivative with respect to φ

def Dfφ(func,φ,δ=1e-4,**kwarg):
    return (func(φ=φ+δ,**kwarg)-func(φ=φ,**kwarg))/δ 


##2nd derivative with respect to φ

def Dfφφ(func,φ,δ=1e-4,**kwarg):
#    value=(func(φ=φ+δ,**kwarg)+func(φ=φ-δ,**kwarg)-2*func(φ=φ,**kwarg))/δ**2 
    return (func(φ=φ+δ,**kwarg)+func(φ=φ-δ,**kwarg)-2*func(φ=φ,**kwarg))/δ**2 
#When we define a function, the arguments can be positional and keyword arguments both at the same time.


#Elctrostatic Interaction Parameters

##Bjerrum length
def lB(T=None):
    return lambda a: a/T

##Debye-Huckel kappa
def kappa(lB=None):
    return lambda *charge_square_ions:np.sqrt(4 * π * lB * sum(charge_square_ions))


#Functions Related to Percolation

##The reacted probability in the sol part
def p_sol(K_c,p_cri=np.inf,tol=1e-7,itn=50,Type="pregel"): #K_c= λφ*f/N
    try:
        p=1+1/(2*K_c)-np.sqrt(np.power(1+1/(2*K_c),2)-1)
        if p > p_cri and Type == "Stockmayer":
            p = p_cri
        elif p > p_cri and Type == "Flory":
            i=1
            f = 1+1/p_cri
            p_s=K_c*np.power(1-p,f) #Input initial vlaue
            while (True):
                i=i+1
                p_ref=p_s
                p_s=K_c*np.power(1-p,f)/np.power(1-p_s,f-2)
                if abs(p_s-p_ref) < tol and i> itn:
                    break
            p=p_s    
        return p
    except:
        return 0
    

##Reaction Rate Constant    
def λ(λ0=None):
    return lambda J, κ=0, g=1: λ0*np.power(exp( J*exp(-κ)*(1+κ)), 2*g-1)
#λ0:collision rate J:association strength g:# of participating monomers κ:inverse DH length    

##Gelation Threshold Concentration
φ_cri=lambda N, λ: 1/λ * (N-1)/np.power(N-2,2)

## single chain monomer concentration: λf/N*φ1=p/(1-p)^{f-2}
φ1=lambda N, φ, λ, gelstate="pregel":( 1/λ * p_sol(λ*φ , 1/(N-1), Type=gelstate) 
                           * np.power(1-p_sol(λ*φ, 1/(N-1),  Type=gelstate), N-2) )

## number concentration of all associations: Σ φ/mN=p/(1-p)^2*(1-pf/2)/λf
d=lambda N, φ, λ, g=1, gelstate="pregel":( p_sol(λ*φ , 1/(N-1), Type=gelstate) 
                             *(1/N - p_sol(λ/g*φ , 1/(N/g-1), Type=gelstate)/2) 
                             /np.power( 1-p_sol(λ*φ,1/(N-1), Type=gelstate),2) / λ )
    

#############################################################    
    
    
class macro_parameters:
    
    
        def __init__(self, **kwarg):   
            # **kwarg includes the frequently used macro-parameters or the invariant parameters during PT 
            self.N  =kwarg.get("N")  #degree of polymerization
            self.T  =kwarg.get("T")  #Temperature
            self.χ  =kwarg.get("χ")  #FH mixiing 
            self.lB =kwarg.get("lB")  #Bjerrum length
            self.α  =kwarg.get("α")  #Degree of ionization
            self.λ0 =kwarg.get("λ0")  #Collision rate
        
        
#############################################################         
    

#Free Energies
#Paradigm
#def free energy(macro-parameters or the invariant parameters during PT):
#    return free energy(order parmeters, tuning parameters)


##Mixing Entropy of Macromolecules
def fmacromix(N=None):
    return lambda φ: φ*log(φ)/N

##Mixing Entropy of molecules
fmolmix=lambda φ0: φ0*log(φ0)    

##Interaction Energy between two speicies. Ex. solvent φ0 and monomer φ
def fχ(χ=None):
    return lambda φ0, φ: χ*φ0*φ

##Free energy due to associations:φ/N*(1+ln(φ1/φ))-Σ φ/mN
def fassociation(N=None):
    return lambda φ, λ, gelstate="pregel": φ/N*(1-log(φ)+log(φ1(N=N, φ=φ, λ=λ, gelstate=gelstate)))-d(N=N, φ=φ, λ=λ, gelstate=gelstate)

##Ion fluctution energy from DH theorem    
fionfluctuation=lambda κ:-np.power(κ,3)/(12*π)

##Long range Kessom dipole-dipole interaction: -4π/18*J^2*φ^2, where J=lB*p^2/l^3 
fdipole_dipole=lambda φ, J, κ=0, ξ=1: ( -π/18*np.power(J,2)/np.power(ξ,3)*np.power(φ,2)
                                                       *exp(-2*κ*ξ)*(4+8*κ*ξ+4*np.power(κ*ξ,2)+np.power(κ*ξ,3)) 
                                                   )
###ξ is the cut-off of dipole forces in the unit of Kuhn length

##Long range Kessom charge-dipole interaction:-π/3 * lB * cion * φ * J
def fcharge_dipole(lB=None): 
    return lambda cion, φ, J, κ, ξ=1: -π/3 * lB * cion * φ * J/ξ * exp(-2*κ*ξ)*(2+κ*ξ)



#The combination way of trees
ω=lambda f,m: np.math.factorial(f*m-m)/(np.math.factorial(m)*np.math.factorial(f*m-2*m+2))
ωxmlarge=lambda f,x,m:( np.sqrt(f-1)/np.power(f-2,5/2)/np.sqrt(2*np.pi)
                  /np.power(m,5/2)*np.power(x/((f-2)**(f-2)/(f-1)**(f-1)),m) ) #Invalid after 3000

#The Sotckmayer's functions S1, S2 ,....
def S(k,f):
    return lambda x:( sum(np.power(m,k)*ω(f,m)*np.power(x,m) for m in range(1, 50))+
                 sum(np.power(m,k)*ωxmlarge(f,x,m) for m in range(51, 500)) )
#The modified Stockmayer's functions
def Sflu(k,f,νf):
    return lambda x,g,h:( 
                        sum(np.power(m,k)*ω(f,m)*np.power(x,m)
                            *np.power(g,m**2-m)/np.power(h,m**(2+2*νf)-m) 
                            for m in range(1, 50))+
                        sum(np.power(m,k)*ωxmlarge(f,x,m)*np.power(g/h,m**2-m)
                            /np.power(h,m**(2+2*νf)-m**(2)) for m in range(51, 500)) )
#weight numbers
def weight(power,**fνf):
    return lambda **xgh: Sflu(k=power+1,**fνf)(**xgh)/Sflu(k=1,**fνf)(**xgh)

#modified factors
def gf(t,Nφ,νfprime,fprime):
    return lambda **xgh: np.exp(np.power(Nφ/weight(power=1+2*νfprime,νf=νfprime,f=fprime)(**xgh),3/2)
                                *np.power(weight(power=1,νf=νfprime,f=fprime)(**xgh),2)
                                *(np.power(t,3*0.63)-0.63*np.power(t,3*0.63-1)/Nφ
                                  /weight(power=1,νf=νfprime,f=fprime)(**xgh))
                                )
    
def hf(t,Nφ,νfprime,fprime):
    return lambda **xgh: np.exp(np.power(Nφ,3/2)*np.power(t,3*0.63)/2
                                *np.power(weight(power=1,νf=νfprime,f=fprime)(**xgh),3)
                                /np.power(weight(power=1+2*νfprime,νf=νfprime,f=fprime)(**xgh),5/2)
                                ) 