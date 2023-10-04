# -*- coding: utf-8 -*-
"""CorrelationFunction.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1Ba0nxdPv8TMGyuT9t-SgNSLXpvtn3aZY
"""

import decimal
from decimal import Decimal
decimal.getcontext().prec = 20

import math
import random
import matplotlib.pyplot as plt
import numpy as np
import sympy
from sympy import symbols, diff
from sympy.combinatorics import Permutation
from sympy.utilities.iterables import multiset_permutations

from numpy.linalg import det,norm

import scipy as scp
from scipy.optimize import minimize,shgo,root
from scipy.misc import derivative
from scipy import integrate
from scipy.linalg import inv,eig

sumk_index=[[0, 1],[0, 1, 0, 2],[0, 1, 0, 2, 0,3]]
precision=1e-8

#Output the position of specific components
def position_kinds(sequence,kinds):       #There are only two kinds in the code    
    return np.where(sequence == kinds)[0]

#The kernals of correlation functions (Ugly in terms of code)  Tolerance=1e-6
def g1(alpha,length,N): #length is the length of bounds
   #alpha=np.float128(alpha)
   alpha=Decimal(alpha)
   length=Decimal(length)
   N=Decimal(N)
   if abs(alpha) < precision:
      value=length/N
   else:
      value=(np.exp(alpha*length)-1)/(N*alpha)
   return value

def g2(alpha,beta,length,N):
    alpha=Decimal(alpha)
    beta=Decimal(beta)
    length=Decimal(length)
    N=Decimal(N)
    if abs(alpha) < precision and abs(beta) < precision:
       value=pow(length/N,2)/2 
    elif abs(alpha) >= precision and abs(beta) < precision:
       value=(length*np.exp(alpha*length))/(pow(N,2)*alpha)-g1(alpha,length,N)/(N*alpha)  
    else:
       value=(g1(alpha+beta,length,N)- g1(alpha,length,N))/(N*beta)
    return value

def g3(alpha,beta,gamma,length,N):
    alpha=Decimal(alpha)
    beta=Decimal(beta)
    gamma=Decimal(gamma)
    length=Decimal(length)
    N=Decimal(N)
    if abs(alpha) < precision and abs(beta) < precision  and  abs(gamma) < precision:
       value= pow(length/N,3)/6
    elif abs(alpha) >= precision and abs(beta) < precision and abs(gamma) < precision:
       value= (pow(length,2)*np.exp(alpha*length))/(2*pow(N,3)*alpha)-g2(alpha,0,length,N)/ (N*alpha)    
    elif abs(beta) >= precision and abs(gamma) < precision:
       value=(g2(alpha+beta,0,length,N)-g2(alpha,beta,length,N))/(N*beta)
    else:      
       value=(g2(alpha,beta+gamma,length,N)- g2(alpha,beta,length,N))/(N*gamma)
    return value

def g4(alpha,beta,gamma,delta,length,N):
    alpha=Decimal(alpha)
    beta=Decimal(beta)
    gamma=Decimal(gamma)
    delta=Decimal(delta)
    length=Decimal(length)
    N=Decimal(N)
    if abs(alpha) < precision and abs(beta) < precision and  abs(gamma) < precision and abs(delta) < precision:
       value=pow(length/N,4)/24
    elif abs(alpha) >= precision and abs(beta) < precision and  abs(gamma) < precision and abs(delta) < precision:
       value=(pow(length,3)*np.exp(alpha*length))/(6*pow(N,4)*alpha)-g3(alpha,0,0,length,N)/ (N*alpha)
    elif  abs(beta) >= precision and abs(gamma) < precision and abs(delta) < precision:
       value=(g3(alpha+beta,0,0,length,N)-g3(alpha,beta,0,length,N))/(N*beta)    
    elif abs(gamma) >= precision and abs(delta) < precision:
       value=(g3(alpha,beta+gamma,0,length,N)-g3(alpha,beta,gamma,length,N))/(N*gamma)
    else:
       value=( g3(alpha,beta,gamma+delta,length,N)- g3(alpha,beta,gamma,length,N) ) / (N*delta)
    return value

#Rank 2 Correlation function
def G2(sequence,kinds,klist,Bounds,N):
  klist=np.array(klist)
  kind1,kind2=kinds
  #Combinations of specific kind domains 
  domainP=np.array(np.meshgrid(position_kinds(sequence,kind1),position_kinds(sequence,kind2))).T.reshape(-1,2) # Incrediblly Simple https://tinyurl.com/2axc5n26
  A1=0
  A2=0

  #Sorting domains according to their positions
  for c in domainP:
      unique,counts=np.unique(c, return_counts=True)  
      if len(counts)== 1:  #2
          per=np.argsort(c)[::-1]
          Bounds_sort=Bounds[c[per]]  #D_(1)= D_(2).......
          D=Decimal(Bounds_sort[0][1]-Bounds_sort[0][0])
          L=Decimal(Bounds_sort[0][0])
          for p in multiset_permutations([0,1]):
              sumk=np.add.reduceat(klist[p],sumk_index[0])[::2]  #k_(1),k_(1)+k_(2).... 
              sumk_square=np.einsum('ij, ij->i', sumk, sumk)
              alpha_i=Decimal(sumk_square[0]/6)
              A1=A1+g2(-alpha_i,alpha_i,D,N)
      else: #1+1
          per=np.argsort(c)[::-1]
          Bounds_sort=Bounds[c[per]]  #D_(1)>= D_(2).......
          sumk=np.add.reduceat(klist[per],sumk_index[0])[::2]  #k_(1),k_(1)+k_(2).... 
          sumk_square=np.einsum('ij, ij->i', sumk, sumk)
          #####################################
          D1=Decimal(Bounds_sort[0][1]-Bounds_sort[0][0])
          L1=Decimal(Bounds_sort[0][0])
          D2=Decimal(Bounds_sort[1][1]-Bounds_sort[1][0])
          L2=Decimal(Bounds_sort[1][0])
          alpha_1=Decimal(sumk_square[0]/6)
          C=Decimal(alpha_1*(L1-L2))
          A2=(A2+np.exp(-C)*g1(-alpha_1,D1,N)*g1(alpha_1,D2,N))
  value=A1+A2
  
  #####
  klist_square=np.einsum('ij, ij->i', klist, klist)
  # if (klist_square<precision).any():
  #     print('Warning:momentums are too small')
  if value<0:
      value='error'
  return value

#Rank 3 Correlation function
def G3(sequence,kinds,klist,Bounds,N):
  klist=np.array(klist)
  kind1,kind2,kind3=kinds
  #Combinations of specific kind domains 
  domainP=np.array(np.meshgrid(
      position_kinds(sequence,kind1),position_kinds(sequence,kind2),position_kinds(sequence,kind3)
      )).T.reshape(-1,3) # Incrediblly Simple https://tinyurl.com/2axc5n26
  A1=0
  A21=0
  A12=0
  A3=0
      
  #Sorting domains according to their positions
  for c in domainP:
      unique,counts=np.unique(c, return_counts=True)  
      if len(counts)==1: #3
          per=np.argsort(c)[::-1]
          Bounds_sort=Bounds[c[per]]
          D=Bounds_sort[0][1]-Bounds_sort[0][0]
          L=Bounds_sort[0][0]
          for p in multiset_permutations([0,1,2]):
            sumk=np.add.reduceat(klist[p],sumk_index[1])[::2]  #k_(1),k_(1)+k_(2).... 
            sumk_square=np.einsum('ij, ij->i', sumk, sumk)
            alpha_i=sumk_square[0]/6
            alpha_ij=sumk_square[1]/6
            A1=A1+g3(-alpha_i,alpha_i-alpha_ij,alpha_ij,D,N)
      elif len(counts)== 2:    
        if (counts == [1, 2]).all():  #smallvalue:1 largevalue:2
          per=np.argsort(c)[::-1]
          Bounds_sort=Bounds[c[per]]
          tempk=klist[per]
          D=Bounds_sort[0][1]-Bounds_sort[0][0]
          L=Bounds_sort[0][0]
          D3=Bounds_sort[2][1]-Bounds_sort[2][0]
          L3=Bounds_sort[2][0]
          for p in multiset_permutations([0,1]):
            sumk=np.add.reduceat(tempk[np.append(p,2)],sumk_index[1])[::2]  #k_(1),k_(1)+k_(2).... 
            sumk_square=np.einsum('ij, ij->i', sumk, sumk)
            alpha_i=sumk_square[0]/6
            alpha_12=sumk_square[1]/6
            C=Decimal(alpha_12*(L-L3))
            A12=A12+np.exp(-C)*g2(-alpha_i,alpha_i-alpha_12,D,N)*g1(alpha_12,D3,N)
        else:    #smallvalue:2 largevalue:1 
          per=np.argsort(c)[::-1]
          Bounds_sort=Bounds[c[per]]
          tempk=klist[per]
          D1=Bounds_sort[0][1]-Bounds_sort[0][0]
          L1=Bounds_sort[0][0]
          D=Bounds_sort[1][1]-Bounds_sort[1][0]
          L=Bounds_sort[1][0]
          for p in multiset_permutations([1,2]):
            sumk=np.add.reduceat(tempk[np.append(0,p)],sumk_index[1])[::2]  #k_(1),k_(1)+k_(2).... 
            sumk_square=np.einsum('ij, ij->i', sumk, sumk)
            alpha_1=sumk_square[0]/6
            alpha_1i=sumk_square[1]/6
            C=Decimal(alpha_1*(L1-L))
            A21=A21+np.exp(-C)*g1(-alpha_1,D1,N)*g2(alpha_1-alpha_1i,alpha_1i,D,N)     
      else: #111
        per=np.argsort(c)[::-1]
        Bounds_sort=Bounds[c[per]]
        D1=Bounds_sort[0][1]-Bounds_sort[0][0]
        L1=Bounds_sort[0][0]
        D2=Bounds_sort[1][1]-Bounds_sort[1][0]
        L2=Bounds_sort[1][0]
        D3=Bounds_sort[2][1]-Bounds_sort[2][0]
        L3=Bounds_sort[2][0]
        sumk=np.add.reduceat(klist[per],sumk_index[1])[::2]  #k_(1),k_(1)+k_(2).... 
        sumk_square=np.einsum('ij, ij->i', sumk, sumk)
        alpha_1=sumk_square[0]/6
        alpha_12=sumk_square[1]/6
        C=Decimal(alpha_1*(L1-L2)+alpha_12*(L2-L3))
        A3=A3+np.exp(-C)*g1(-alpha_1,D1,N)*g1(alpha_1-alpha_12,D2,N)*g1(alpha_12,D3,N)
  value=A1+A12+A21+A3


  ####
  klist_square=np.einsum('ij, ij->i', klist, klist)
  # if (klist_square<precision).any():
  #     print('Warning:momentums are too small')
  if value<0:
      value='error'
  return value

#Rank 4 Correlation function
def G4(sequence,kinds,klist,Bounds,N):
  klist=np.array(klist)
  kind1,kind2,kind3,kind4=kinds
  #Combinations of specific kind domains 
  domainP=np.array(np.meshgrid(
      position_kinds(sequence,kind1),position_kinds(sequence,kind2),position_kinds(sequence,kind3),position_kinds(sequence,kind4)
      )).T.reshape(-1,4) # Incrediblly Simple https://tinyurl.com/2axc5n26
  A1=0
  A13=0
  A31=0
  A22=0
  A112=0
  A121=0
  A211=0
  A4=0
      
  #Sorting domains according to their positions
  for c in domainP:
      unique,counts=np.unique(c, return_counts=True)  
      if len(counts)== 1:
        per=np.argsort(c)[::-1]
        Bounds_sort=Bounds[c[per]]
        D=Bounds_sort[0][1]-Bounds_sort[0][0]
        L=Bounds_sort[0][0]
        for p in multiset_permutations([0,1,2,3]):
            sumk=np.add.reduceat(klist[p],sumk_index[2])[::2]  #k_(1),k_(1)+k_(2).... 
            sumk_square=np.einsum('ij, ij->i', sumk, sumk)
            alpha_i=sumk_square[0]/6
            alpha_ij=sumk_square[1]/6
            alpha_ijk=sumk_square[2]/6
            A1=A1+g4(-alpha_i,alpha_i-alpha_ij,alpha_ij-alpha_ijk,alpha_ijk,D,N)
      elif len(counts)== 2:
          if (counts == [2, 2]).all(): #smallvalue:2 largevalue:2
            per=np.argsort(c)[::-1]
            Bounds_sort=Bounds[c[per]]
            Dprime=Bounds_sort[0][1]-Bounds_sort[0][0]
            Lprime=Bounds_sort[0][0]
            D=Bounds_sort[2][1]-Bounds_sort[2][0]
            L=Bounds_sort[2][0]
            tempk=klist[per]
            for p in multiset_permutations([0,1]):
              for q in multiset_permutations([2,3]):
                sumk=np.add.reduceat(tempk[np.append(p,q)],sumk_index[2])[::2]  #k_(1),k_(1)+k_(2).... 
                sumk_square=np.einsum('ij, ij->i', sumk, sumk)
                alpha_i=sumk_square[0]/6
                alpha_12=sumk_square[1]/6
                alpha_12j=sumk_square[2]/6
                C=Decimal(alpha_12*(Lprime-L))
                A22=A22+np.exp(-C)*g2(-alpha_i,alpha_i-alpha_12,Dprime,N)*g2(alpha_12-alpha_12j,alpha_12j,D,N)
          elif (counts == [1, 3]).all(): #smallvalue:1 largevalue:3
            per=np.argsort(c)[::-1]
            Bounds_sort=Bounds[c[per]]
            D=Bounds_sort[0][1]-Bounds_sort[0][0]
            L=Bounds_sort[0][0]
            D4=Bounds_sort[3][1]-Bounds_sort[3][0]
            L4=Bounds_sort[3][0]
            tempk=klist[per]
            for p in multiset_permutations([0,1,2]):
              sumk=np.add.reduceat(tempk[np.append(p,3)],sumk_index[2])[::2]  #k_(1),k_(1)+k_(2).... 
              sumk_square=np.einsum('ij, ij->i', sumk, sumk)
              alpha_i=sumk_square[0]/6
              alpha_ij=sumk_square[1]/6
              alpha_123=sumk_square[2]/6
              C=Decimal(alpha_123*(L-L4))
              A13=A13+np.exp(-C)*g3(-alpha_i,alpha_i-alpha_ij,alpha_ij-alpha_123,D,N)*g1(alpha_123,D4,N)
          else: #smallvalue:3 largevalue:1
            per=np.argsort(c)[::-1]
            Bounds_sort=Bounds[c[per]]
            D1=Bounds_sort[0][1]-Bounds_sort[0][0]
            L1=Bounds_sort[0][0]
            D=Bounds_sort[1][1]-Bounds_sort[1][0]
            L=Bounds_sort[1][0]
            tempk=klist[per]
            for p in multiset_permutations([1,2,3]):
              sumk=np.add.reduceat(tempk[np.append(0,p)],sumk_index[2])[::2]  #k_(1),k_(1)+k_(2).... 
              sumk_square=np.einsum('ij, ij->i', sumk, sumk)
              alpha_1=sumk_square[0]/6
              alpha_1i=sumk_square[1]/6
              alpha_1ij=sumk_square[2]/6
              C=Decimal(alpha_1*(L1-L))
              A31=A31+np.exp(-C)*g1(-alpha_1,D1,N)*g3(alpha_1-alpha_1i,alpha_1i-alpha_1ij,alpha_1ij,D,N)
      elif len(counts)==3:
          if (counts == [1, 1, 2]).all(): #small:1 medium:1 large:2
            per=np.argsort(c)[::-1]
            Bounds_sort=Bounds[c[per]]
            D=Bounds_sort[0][1]-Bounds_sort[0][0]
            L=Bounds_sort[0][0]
            D3=Bounds_sort[2][1]-Bounds_sort[2][0]
            L3=Bounds_sort[2][0]
            D4=Bounds_sort[3][1]-Bounds_sort[3][0]
            L4=Bounds_sort[3][0]
            tempk=klist[per]
            for p in multiset_permutations([0,1]):
              sumk=np.add.reduceat(tempk[np.append(p,[2,3])],sumk_index[2])[::2]  #k_(1),k_(1)+k_(2).... 
              sumk_square=np.einsum('ij, ij->i', sumk, sumk)
              alpha_i=sumk_square[0]/6
              alpha_12=sumk_square[1]/6
              alpha_123=sumk_square[2]/6
              C=Decimal(alpha_12*(L-L3)+alpha_123*(L3-L4))
              A112=A112+np.exp(-C)*g2(-alpha_i,alpha_i-alpha_12,D,N)*g1(alpha_12-alpha_123,D3,N)*g1(alpha_123,D4,N)
          elif (counts == [1, 2, 1]).all():#small:1 medium:2 large:1
            per=np.argsort(c)[::-1]
            Bounds_sort=Bounds[c[per]]
            D1=Bounds_sort[0][1]-Bounds_sort[0][0]
            L1=Bounds_sort[0][0]
            D=Bounds_sort[1][1]-Bounds_sort[1][0]
            L=Bounds_sort[1][0]
            D4=Bounds_sort[3][1]-Bounds_sort[3][0]
            L4=Bounds_sort[3][0]
            tempk=klist[per]
            for p in multiset_permutations([1,2]):
              sumk=np.add.reduceat(tempk[np.append(np.append(0,p),3)],sumk_index[2])[::2]  #k_(1),k_(1)+k_(2).... 
              sumk_square=np.einsum('ij, ij->i', sumk, sumk)
              alpha_1=sumk_square[0]/6
              alpha_1i=sumk_square[1]/6
              alpha_123=sumk_square[2]/6
              C=Decimal(alpha_1*(L1-L)+alpha_123*(L-L4))
              A121=A121+np.exp(-C)*g1(-alpha_1,D1,N)*g2(alpha_1-alpha_1i,alpha_1i-alpha_123,D,N)*g1(alpha_123,D4,N)
          else: #small:2 medium:1 large:1
            per=np.argsort(c)[::-1]
            Bounds_sort=Bounds[c[per]]
            D1=Bounds_sort[0][1]-Bounds_sort[0][0]
            L1=Bounds_sort[0][0]
            D2=Bounds_sort[1][1]-Bounds_sort[1][0]
            L2=Bounds_sort[1][0]
            D=Bounds_sort[2][1]-Bounds_sort[2][0]
            L=Bounds_sort[2][0]
            tempk=klist[per]
            for p in multiset_permutations([2,3]):
              sumk=np.add.reduceat(tempk[np.append([0,1],p)],sumk_index[2])[::2]  #k_(1),k_(1)+k_(2).... 
              sumk_square=np.einsum('ij, ij->i', sumk, sumk)
              alpha_1=sumk_square[0]/6
              alpha_12=sumk_square[1]/6
              alpha_12i=sumk_square[2]/6
              C=Decimal(alpha_1*(L1-L2)+alpha_12*(L2-L))
              A211=A211+np.exp(-C)*g1(-alpha_1,D1,N)*g1(alpha_1-alpha_12,D2,N)*g2(alpha_12-alpha_12i,alpha_12i,D,N)
      else:
        per=np.argsort(c)[::-1]
        Bounds_sort=Bounds[c[per]]
        D1=Bounds_sort[0][1]-Bounds_sort[0][0]
        L1=Bounds_sort[0][0]
        D2=Bounds_sort[1][1]-Bounds_sort[1][0]
        L2=Bounds_sort[1][0]
        D3=Bounds_sort[2][1]-Bounds_sort[2][0]
        L3=Bounds_sort[2][0]
        D4=Bounds_sort[3][1]-Bounds_sort[3][0]
        L4=Bounds_sort[3][0]
        sumk=np.add.reduceat(klist[per],sumk_index[2])[::2]  #k_(1),k_(1)+k_(2).... 
        sumk_square=np.einsum('ij, ij->i', sumk, sumk)
        alpha_1=sumk_square[0]/6
        alpha_12=sumk_square[1]/6
        alpha_123=sumk_square[2]/6
        C=Decimal(alpha_1*(L1-L2)+alpha_12*(L2-L3)+alpha_123*(L3-L4))
        A4=A4+np.exp(-C)*g1(-alpha_1,D1,N)*g1(alpha_1-alpha_12,D2,N)*g1(alpha_12-alpha_123,D3,N)*g1(alpha_123,D4,N)
  value=A1+A22+A13+A31+A112+A121+A211+A4
  ####
  klist_square=np.einsum('ij, ij->i', klist, klist)
  # if (klist_square<precision).any():
  #     print('Warning:momentums are too small')
  if value<0:
     value='error'
  return value

"""Vertex function of polymer"""

#!!only work for 2-kinds monomer
dim=2
def Gamma2_p(seq,phi_p,q,Bounds,N):
  if q !=0 :
    k1=np.array([q,0,0])
    k2=np.array([-q,0,0])
    klist=np.array([k1,k2])
    S2_p=np.zeros((dim,dim))
    for i in range(0,dim):
     for j in range(0,dim):
      S2_p[i,j]=N*phi_p*float(G2(seq,[i,j],klist,Bounds,N))
    matrix = inv(S2_p)
  else:
    matrix=np.ones((dim,dim))*1/(N*phi_p)
    #np.array([[1/(N*phi_p),1/(N*phi_p)],[1/(N*phi_p),1/(N*phi_p)]])  
  return matrix

def Gamma3_p(seq,phi_p,klist,Bounds,N):
    q1,q2,q3=klist
    S3_p=np.zeros((dim,dim,dim))
    for i in range(0,dim):
      for j in range(0,dim):
        for k in range(0,dim):
           S3_p[i,j,k]=N**2*phi_p*float(G3(seq,[i,j,k],klist,Bounds,N))             
    matrix=np.einsum('ix,jy,kz,xyz->ijk',Gamma2_p(seq,phi_p,norm(q1),Bounds,N),Gamma2_p(seq,phi_p,norm(q2),Bounds,N),Gamma2_p(seq,phi_p,norm(q3),Bounds,N),S3_p)
    return matrix

def Gamma4_p(seq,phi_p,klist,Bounds,N):
    q1,q2,q3,q4=klist
    matrix=np.zeros((dim,dim,dim,dim))
    S4_p=np.zeros((dim,dim,dim,dim))
    S4_p1p2=np.zeros((dim,dim,dim,dim)) #Gamma3 G2 Gamma3
    S4_p1p3=np.zeros((dim,dim,dim,dim))
    S4_p1p4=np.zeros((dim,dim,dim,dim))
    for i in range(0,dim):
      for j in range(0,dim):
        for k in range(0,dim):
          for l in range(0,dim):
           S4_p[i,j,k,l]=N**3*phi_p*float(G4(seq,[i,j,k,l],klist,Bounds,N))
           for x in range(0,dim):
             for y in range(0,dim):
                  S4_p1p2[i,j,k,l]=(S4_p1p2[i,j,k,l]+
                  (N**2*phi_p*float(G3(seq,[i,j,x],[q1,q2,-q1-q2],Bounds,N))
                  *Gamma2_p(seq,phi_p,norm(q1+q2),Bounds,N)[x,y]
                  *N**2*phi_p*float(G3(seq,[y,k,l],[q1+q2,q3,q4],Bounds,N))))

                  S4_p1p3[i,j,k,l]=(S4_p1p3[i,j,k,l]+
                  (N**2*phi_p*float(G3(seq,[i,k,x],[q1,q3,-q1-q3],Bounds,N))
                  *Gamma2_p(seq,phi_p,norm(q1+q3),Bounds,N)[x,y]
                  *N**2*phi_p*float(G3(seq,[y,j,l],[q1+q3,q2,q4],Bounds,N))))

                  S4_p1p4[i,j,k,l]=(S4_p1p4[i,j,k,l]+
                  (N**2*phi_p*float(G3(seq,[i,l,x],[q1,q4,-q1-q4],Bounds,N))
                  *Gamma2_p(seq,phi_p,norm(q1+q4),Bounds,N)[x,y]
                  *N**2*phi_p*float(G3(seq,[y,k,j],[q1+q4,q3,q2],Bounds,N))))
    #print(S4_p1p2+S4_p1p3+S4_p1p4)
    matrix_1=S4_p1p2+S4_p1p3+S4_p1p4 - S4_p                
    matrix=np.einsum('iw,jx,ky,lz,wxyz->ijkl',Gamma2_p(seq,phi_p,norm(q1),Bounds,N),Gamma2_p(seq,phi_p,norm(q2),Bounds,N),
                       Gamma2_p(seq,phi_p,norm(q3),Bounds,N),Gamma2_p(seq,phi_p,norm(q4),Bounds,N),matrix_1)
    return matrix