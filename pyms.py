# -*- coding: utf-8 -*-
import sys
import pdb
import math
import pymsfunc as fn
import numpy
import cmath
import scipy
from scipy import linalg
import matplotlib 
from matplotlib import pylab as pt
import matplotlib.pyplot as pt
from pylab import *
import time


def main():
    
    solar_point=open('sollo20.dat','r')
    solar_dat=solar_point.readlines()
    
    mass_point=open('isotope_As.dat','r')
    isomass=mass_point.readlines()
    
    ion=[]*len(solar_dat)
    abu_array=[]*len(solar_dat)
    for x in range(len(solar_dat)):
      s=solar_dat[x].split()
      ion.append(s[0])
      abu_array.append(s[1])
    
    for x in range(len(abu_array)):
      abu_array[x]=float(abu_array[x])
      isomass[x]=float(isomass[x])

    sion=fn.main_s_target()
    dion=fn.main_s_daughter()
    uniq=fn.uniqmake(sion)
   
    sig=fn.sigs(sion)  
    
    b=fn.main_branch(sion)
    mass=fn.massmake()
   
    Se79=fn.where(sion,"Se79")
    Zr93=fn.where(uniq,"Zr93")
    Tc99=fn.where(uniq,"Tc99")
    Pd107=fn.where(uniq,"Pd107")
    I129=fn.where(uniq,"I129")
    Pb205=fn.where(uniq,"Pb205")

    Br79=fn.where(ion,"Br79")
    Nb93=fn.where(ion,"Nb93")
    Ru99=fn.where(ion,"Ru99")
    Ag107=fn.where(ion,"Ag107")
    Xe129=fn.where(ion,"Xe129")
    Tl205=fn.where(ion,"Tl205")
    
    uniq2=copy(uniq)
    uniq2[Se79]="Br79"
    uniq2[Zr93]="Nb93"
    uniq2[Tc99]="Ru99"
    uniq2[Pd107]="Ag107"
    uniq2[I129]="Xe129"
    uniq2[Pb205]="Tl205"
    
    s_index=fn.index(uniq2,ion)
    
   
    solar=[0]*len(uniq2)
    for x in range(len(s_index)):
      solar[x]=abu_array[s_index[x]]
    Bi209=fn.where(ion,"Bi209")
    solar[len(uniq2)-1]=abu_array[Bi209]
   
    '''
    solar=[0]*len(uniq)
    for x in range(len(uniq)):
      if x==Se79: solar[x]=abu_array[Br79]
      elif x==Zr93: solar[x]=abu_array[Nb93] 
      elif x==Tc99: solar[x]=abu_array[Ru99]
      elif x==Pd107: solar[x]=abu_array[Ag107]
      elif x==I129: solar[x]=abu_array[Xe129]
      elif x==Pb205: solar[x]=abu_array[Tl205]
      else: solar[x]=abu_array[x]
      
  
     
    #sig=numpy.ones(len(sion),dtype=np.float)
    #BEGIN TAU LOOP
    
    upp=50
    
    uniq=uniq[0:upp]
    sion=sion[0:upp]
    dion=dion[0:upp]
    sig=sig[0:upp]
    b=b[0:upp]
    solar=solar[0:upp]
    mass=mass[0:upp]
    '''
    
    weak=0
    if weak==1:
      sion=sion[0:48]
      dion=dion[0:48]
      sig=sig[0:48]
      b=b[0:48]
      uniq=uniq[0:48]
      solar=solar[0:48]
      mass=mass[0:48]
      uniq2=uniq2[0:48]
    
    #print sion
    #print uniq
    num=float(20)
    up=float(0.5)
    down=float(0.01)
    tau=fn.makearray(num,up,down)
    #tau[0.277]
    tau=[1.24,.3]
    aa=zeros([4,len(tau)],float)
    bb=ones([4,1],float)
    sfits=zeros([len(solar),len(tau)],float)
    for k in range(len(tau)):
      
      a=fn.matrix(uniq,sion,dion,sig,b)   
      #a=(zip(*a)) 
      #print a[40:49]
      a=numpy.multiply(a,tau[k])
      aexp=scipy.linalg.expm(a)
    
      n=[0]*len(uniq)
      for x in range(len(n)):
        n[x]=aexp[x][0] 
      

      for x in range(len(n)):
        if n[x] <= 0.0: 
          print(uniq[x],x)
      #convert n to mass frac
      n=fn.atoms2massfrac(n,solar,abu_array)
      #for x in range(len(n)):
	#if n[x] > 0.0:n[x]=math.log10(n[x]/solar[x])
      #pt.plot(mass,n,'*')
      #pt.show()
      #print aexp
     
      ratio=[0]*len(n)
      for x in range(len(n)):
	      ratio[x]=n[x]/solar[x]
      
      fac=max(ratio)
      maxx=fn.where(ratio,fac)
      print(uniq2[maxx],tau[k])
      
      for x in range(len(n)):
	      n[x]=n[x]/fac
	      if n[x] <= 0.0:
          print(uniq2[x])
        
        for x in range(len(n)):
  #if ratio[x]> 0.0:ratio[x]=math.log10(n[x]/solar[x])
	        ratio[x]=(n[x]/solar[x])
	#elif ratio[x]< 0.0:print uniq[x]
	        if n[x] > 0.0:
            n[x]=math.log10(n[x])
      
      tit=uniq2[maxx]+' '+str(tau[k])
      print tit
      pt.plot(mass,ratio,'*')
      pt.title(tit)
      pt.ylabel('Log(sproc/solar)')
      pt.xlabel('Mass Number')
      #pt.show()
      
      print tau[k]
      pt.close('all')
      
      Ge70=fn.where(uniq,'Ge70')
      Se76=fn.where(uniq,'Se76')
      Kr82=fn.where(uniq,'Kr82')
      Sr86=fn.where(uniq,'Sr86')
      
      aa[0,k]=ratio[Ge70]
      aa[1,k]=ratio[Se76]
      aa[2,k]=ratio[Kr82]
      aa[3,k]=ratio[Sr86]
      
      for x in range(len(solar)):
	      sfits[x,k]=ratio[x]
        
    #END TAU LOOP
    
    array=zeros(len(solar))
    cc=numpy.linalg.lstsq(aa,bb)
    for i in range(len(tau)):
      for j in range(len(solar)):
	      array[j]+=cc[0][i]*sfits[j,i]

    for x in range(len(solar)):
      if array[x]>=0.0:array[x]=math.log10(array[x])
    
    print(cc[0])
    pt.plot(mass,array,'*')
    pt.show()
    
    # expm --- matrix exponential using Pade approx.
    # expm2 --- matrix exponential using Eigenvalue decomp.
    # expm3 --- matrix exponential using Taylor-series expansion
   
    #aexp=expm(a)
    #n=aexp[*][0]
    
  
main()  