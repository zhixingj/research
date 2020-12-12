# -*- coding: utf-8 -*-

import sys
import pdb
import math

####################################################################
def main_s_target():
  sion=["Fe56",
  "Fe57", 
  "Fe58", 
  "Co59", 
  "Ni60", 
  "Ni61", 
  "Ni62", 
  "Cu63", 
  "Cu63",
  "Ni64",  
  "Zn64", 
  "Cu65",
  "Zn66", 
  "Zn67", 
  "Zn68", 
  "Ga69", 
  "Ge70",
  "Ga71", 
  "Ge72", 
  "Ge73", 
  "Ge74", 
  "As75", 
  "Se76",
  "Se77", 
  "Se78", 
  "Se79", 
  "Se80", 
  "Br81", 
  "Kr82", 
  "Kr83", 
  "Kr84", 
  "Rb85", 
  "Sr86", 
  "Sr87", 
  "Sr88", 
  "Y89", 
  "Zr90", 
  "Zr91", 
  "Zr92", 
  "Zr93", 
  "Zr94", 
  "Mo95", 
  "Mo96", 
  "Mo97", 
  "Mo98", 
  "Tc99", 
  "Ru100", 
  "Ru101",
  "Ru102",
  "Rh103",
  "Pd104",
  "Pd105",
  "Pd106",
  "Pd107",
  "Pd108",
  "Ag109",
  "Cd110",
  "Cd111",
  "Cd112",
  "Cd113",
  "Cd114",
  "In115",
  "Sn116",
  "Sn117",
  "Sn118",
  "Sn119",
  "Sn120",
  "Sb121",
  "Te122",
  "Te123",
  "Te124",
  "Te125",
  "Te126",
  "I127",
  "I127",
  "Te128",
  "I129",
  "Xe128",
  "Xe129",
  "Xe130",
  "Xe131",
  "Xe132",
  "Cs133",
  "Ba134",
  "Ba135",
  "Ba136",
  "Ba137",
  "Ba138",
  "La139",
  "Ce140",
  "Pr141",
  "Nd142",
  "Nd143",
  "Nd144",
  "Nd145",
  "Nd146",
  "Sm147",
  "Sm148",
  "Sm149",
  "Sm150",
  "Eu151",
  "Eu151",
  "Sm152",
  "Gd152",
  "Eu153",
  "Gd154",
  "Gd155",
  "Gd156",
  "Gd157",
  "Gd158",
  "Tb159",
  "Dy160",
  "Dy161",
  "Dy162",
  "Dy163",
  "Dy164",
  "Ho165",
  "Er166",
  "Er167",
  "Er168",
  "Tm169",
  "Yb170",
  "Yb171",
  "Yb172",
  "Yb173",
  "Yb174",
  "Lu175",
  "Lu176",
  "Hf177",
  "Hf178",
  "Hf179",
  "Hf180",
  "Ta181",
  "W182",
  "W183",
  "W184",
  "Re185",
  "Re185",
  "W186",
  "Re187",
  "Os186",
  "Os187",
  "Os188",
  "Os189",
  "Os190",
  "Ir191",
  "Ir191",
  "Os192",
  "Pt192",
  "Ir193",
  "Pt194",
  "Pt195",
  "Pt196",
  "Au197",
  "Hg198",
  "Hg199",
  "Hg200",
  "Hg201",
  "Hg202",
  "Tl203",
  "Pb204",
  "Pb205",
  "Pb206",
  "Pb207",
  "Pb208",
  "Bi209"]
  return sion
####################################################################
def main_s_daughter():
  dion=["Fe57", 
  "Fe58", 
  "Co59", 
  "Ni60", 
  "Ni61", 
  "Ni62", 
  "Cu63", 
  "Ni64",  
  "Zn64", 
  "Cu65",
  "Cu65",
  "Zn66", 
  "Zn67", 
  "Zn68", 
  "Ga69", 
  "Ge70",
  "Ga71", 
  "Ge72", 
  "Ge73", 
  "Ge74", 
  "As75", 
  "Se76",
  "Se77", 
  "Se78", 
  "Se79", 
  "Se80", 
  "Br81", 
  "Kr82", 
  "Kr83", 
  "Kr84", 
  "Rb85", 
  "Sr86", 
  "Sr87", 
  "Sr88", 
  "Y89", 
  "Zr90", 
  "Zr91", 
  "Zr92", 
  "Zr93", 
  "Zr94", 
  "Mo95", 
  "Mo96", 
  "Mo97", 
  "Mo98", 
  "Tc99", 
  "Ru100", 
  "Ru101",
  "Ru102",
  "Rh103",
  "Pd104",
  "Pd105",
  "Pd106",
  "Pd107",
  "Pd108",
  "Ag109",
  "Cd110",
  "Cd111",
  "Cd112",
  "Cd113",
  "Cd114",
  "In115",
  "Sn116",
  "Sn117",
  "Sn118",
  "Sn119",
  "Sn120",
  "Sb121",
  "Te122",
  "Te123",
  "Te124",
  "Te125",
  "Te126",
  "I127",
  "Te128",
  "Xe128",
  "I129",
  "Xe130",
  "Xe129",
  "Xe130",
  "Xe131",
  "Xe132",
  "Cs133",
  "Ba134",
  "Ba135",
  "Ba136",
  "Ba137",
  "Ba138",
  "La139",
  "Ce140",
  "Pr141",
  "Nd142",
  "Nd143",
  "Nd144",
  "Nd145",
  "Nd146",
  "Sm147",
  "Sm148",
  "Sm149",
  "Sm150",
  "Eu151",
  "Sm152",
  "Gd152",
  "Eu153",
  "Eu153",
  "Gd154",
  "Gd155",
  "Gd156",
  "Gd157",
  "Gd158",
  "Tb159",
  "Dy160",
  "Dy161",
  "Dy162",
  "Dy163",
  "Dy164",
  "Ho165",
  "Er166",
  "Er167",
  "Er168",
  "Tm169",
  "Yb170",
  "Yb171",
  "Yb172",
  "Yb173",
  "Yb174",
  "Lu175",
  "Lu176",
  "Hf177",
  "Hf178",
  "Hf179",
  "Hf180",
  "Ta181",
  "W182",
  "W183",
  "W184",
  "Re185",
  "W186",
  "Os186",
  "Re187",
  "Os188",
  "Os187",
  "Os188",
  "Os189",
  "Os190",
  "Ir191",
  "Os192",
  "Pt192",
  "Ir193",
  "Ir193",
  "Pt194",
  "Pt195",
  "Pt196",
  "Au197",
  "Hg198",
  "Hg199",
  "Hg200",
  "Hg201",
  "Hg202",
  "Tl203",
  "Pb204",
  "Pb205",
  "Pb206",
  "Pb207",
  "Pb208",
  "Bi209",
  "Pb206"]
  return dion
##################################################################################
def uniqmake(seq, idfun=None):
  result=["Fe56",
  "Fe57", 
  "Fe58", 
  "Co59", 
  "Ni60", 
  "Ni61", 
  "Ni62", 
  "Cu63", 
  "Ni64",  
  "Zn64", 
  "Cu65",
  "Zn66", 
  "Zn67", 
  "Zn68", 
  "Ga69", 
  "Ge70",
  "Ga71", 
  "Ge72", 
  "Ge73", 
  "Ge74", 
  "As75", 
  "Se76",
  "Se77", 
  "Se78", 
  "Se79", 
  "Se80", 
  "Br81", 
  "Kr82", 
  "Kr83", 
  "Kr84", 
  "Rb85", 
  "Sr86", 
  "Sr87", 
  "Sr88", 
  "Y89", 
  "Zr90", 
  "Zr91", 
  "Zr92", 
  "Zr93", 
  "Zr94", 
  "Mo95", 
  "Mo96", 
  "Mo97", 
  "Mo98", 
  "Tc99", 
  "Ru100", 
  "Ru101",
  "Ru102",
  "Rh103",
  "Pd104",
  "Pd105",
  "Pd106",
  "Pd107",
  "Pd108",
  "Ag109",
  "Cd110",
  "Cd111",
  "Cd112",
  "Cd113",
  "Cd114",
  "In115",
  "Sn116",
  "Sn117",
  "Sn118",
  "Sn119",
  "Sn120",
  "Sb121",
  "Te122",
  "Te123",
  "Te124",
  "Te125",
  "Te126",
  "I127",
  "Te128",
  "I129",
  "Xe128",
  "Xe129",
  "Xe130",
  "Xe131",
  "Xe132",
  "Cs133",
  "Ba134",
  "Ba135",
  "Ba136",
  "Ba137",
  "Ba138",
  "La139",
  "Ce140",
  "Pr141",
  "Nd142",
  "Nd143",
  "Nd144",
  "Nd145",
  "Nd146",
  "Sm147",
  "Sm148",
  "Sm149",
  "Sm150",
  "Eu151",
  "Sm152",
  "Gd152",
  "Eu153",
  "Gd154",
  "Gd155",
  "Gd156",
  "Gd157",
  "Gd158",
  "Tb159",
  "Dy160",
  "Dy161",
  "Dy162",
  "Dy163",
  "Dy164",
  "Ho165",
  "Er166",
  "Er167",
  "Er168",
  "Tm169",
  "Yb170",
  "Yb171",
  "Yb172",
  "Yb173",
  "Yb174",
  "Lu175",
  "Lu176",
  "Hf177",
  "Hf178",
  "Hf179",
  "Hf180",
  "Ta181",
  "W182",
  "W183",
  "W184",
  "Re185",
  "W186",
  "Re187",
  "Os186",
  "Os187",
  "Os188",
  "Os189",
  "Os190",
  "Ir191",
  "Os192",
  "Pt192",
  "Ir193",
  "Pt194",
  "Pt195",
  "Pt196",
  "Au197",
  "Hg198",
  "Hg199",
  "Hg200",
  "Hg201",
  "Hg202",
  "Tl203",
  "Pb204",
  "Pb205",
  "Pb206",
  "Pb207",
  "Pb208",
  "Bi209"]
  return result
###############################################################################3
def massmake():
  mass=[56.0,
  57.0,
  58.0,
  59.0,
  60.0,
  61.0,
  62.0,
  63.0,
  64.0,
  64.0,
  65.0,
  66.0,
  67.0,
  68.0,
  69.0,
  70.0,
  71.0,
  72.0,
  73.0,
  74.0,
  75.0,
  76.0,
  77.0,
  78.0,
  79.0,
  80.0,
  81.0,
  82.0,
  83.0,
  84.0,
  85.0,
  86.0,
  87.0,
  88.0,
  89.0,
  90.0,
  91.0,
  92.0,
  93.0,
  94.0,
  95.0,
  96.0,
  97.0,
  98.0,
  99.0,
  100.0,
  101.0,
  102.0,
  103.0,
  104.0,
  105.0,
  106.0,
  107.0,
  108.0,
  109.0,
  110.0,
  111.0,
  112.0,
  113.0,
  114.0,
  115.0,
  116.0,
  117.0,
  118.0,
  119.0,
  120.0,
  121.0,
  122.0,
  123.0,
  124.0,
  125.0,
  126.0,
  127.0,
  128.0,
  129.0,
  128.0,
  129.0,
  130.0,
  131.0,
  132.0,
  133.0,
  134.0,
  135.0,
  136.0,
  137.0,
  138.0,
  139.0,
  140.0,
  141.0,
  142.0,
  143.0,
  144.0,
  145.0,
  146.0,
  147.0,
  148.0,
  149.0,
  150.0,
  151.0,
  152.0,
  152.0,
  153.0,
  154.0,
  155.0,
  156.0,
  157.0,
  158.0,
  159.0,
  160.0,
  161.0,
  162.0,
  163.0,
  164.0,
  165.0,
  166.0,
  167.0,
  168.0,
  169.0,
  170.0,
  171.0,
  172.0,
  173.0,
  174.0,
  175.0,
  176.0,
  177.0,
  178.0,
  179.0,
  180.0,
  181.0,
  182.0,
  183.0,
  184.0,
  185.0,
  186.0,
  186.0,
  187.0,
  187.0,
  188.0,
  189.0,
  190.0,
  191.0,
  192.0,
  192.0,
  193.0,
  194.0,
  195.0,
  196.0,
  197.0,
  198.0,
  199.0,
  200.0,
  201.0,
  202.0,
  203.0,
  204.0,
  205.0,
  206.0,
  207.0,
  208.0,
  209.0]
  return mass
#####################################################################################
def index(uniq,ion):
  s_index=[]*len(uniq)
  for x in range(len(uniq)):
    for y in range(len(ion)):
      if uniq[x]==ion[y]:
	      s_index.append(y)
	
  for x in range(len(s_index)):
      s_index[x]=int(s_index[x])
  return s_index
#####################################################################################
def sigs(sion):
  kad=1
  kep=0
  if kad==1:
    macs_point=open('py_kadonis_macs.dat','r')
    macs_dat=macs_point.readlines()
    
    macs_n=[]*len(macs_dat) #names
    macs_s=[]*len(macs_dat) #values
    for x in range(len(macs_dat)):
      s=macs_dat[x].split()
      macs_n.append(s[0])
      macs_s.append(s[1])
      
    for x in range(len(macs_s)):
      macs_s[x]=float(macs_s[x])
    
    sig=[0]*len(sion)
    for x in range(len(sion)):
      for y in range(len(macs_s)):
        if sion[x]==macs_n[y]:
          sig[x]=macs_s[y]
          return sig
  if kep==1:
    keps_point=open('py_rates.dat','r')
    keps_dat=keps_point.readlines()
    
    keps_n=[]*len(keps_dat)
    keps_s=[]*len(keps_dat)
    
    for x in range(len(keps_dat)):
      s=keps_dat[x].split()
      keps_n.append(s[0])
      keps_s.append(s[1])
      
    for x in range(len(keps_s)):
      keps_s[x]=float(keps_s[x])
    
    sig=[0]*len(sion)
    for x in range(len(sion)):
      for y in range(len(keps_s)):
	if sion[x]==keps_n[y]:
	   sig[x]=keps_s[y]
    
    #normalize rates to Fe56 macs (11.7 mb)
    fac=float(11.7)/sig[0]
    for x in range(len(sig)):
      sig[x]=sig[x]*fac
    return sig
    
#####################################################################################
def main_branch(sion):
  b=[0]*len(sion)
  for x in range(len(sion)):
    b[x]=float(1)
  #branchings (defined by b+)

  Cu63=sion.index("Cu63")
  b[Cu63]=float(0.39)
  b[Cu63+1]=float(1)-b[Cu63]
  
  I127=sion.index("I127")
  b[I127]=float(0.069)
  b[I127+1]=float(1)-b[I127]

  Eu151=sion.index("Eu151")
  b[Eu151]=float(0.2792)
  b[Eu151+1]=float(1)-b[Eu151]

  Re185=sion.index("Re185")
  b[Re185]=float(0.069)
  b[Re185+1]=float(1)-b[Re185]

  Ir191=sion.index("Ir191")
  b[Ir191]=float(0.0476)
  b[Ir191+1]=float(1)-b[Ir191]
  return b
#####################################################################################
def unique(seq, idfun=None):
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        # in old Python versions:
        # if seen.has_key(marker)
        # but in new ones:
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result
#####################################################################################    
def matrix(uniq,sion,dion,sig,b):
  #CONSTRUCTS COEFFICIENT MATRIX
  #ASSUMES THE USE OF COLUMN VECTORS FOR dN/dt = A*N
  #IF DESIRE ROW VECTORS (ARRAYS) TAKE TRANSPOSE

  a=[[0 for col in range(len(uniq))] for row in range(len(uniq))]
  
  for x in range(len(sion)):
    for y in range(len(uniq)):
      if sion[x]==uniq[y]:
	a[y][y]=float(-1)*sig[x]*b[x] + float(-1)*sig[x]*(float(1)-b[x])
	
  for x in range(len(dion)):
    for y in range(len(uniq)):
      if dion[x]==uniq[y]:
	c=y
    for z in range(len(uniq)):
      if sion[x]==uniq[z]:
	r=z
    if c!=r:a[c][r]=sig[x]*b[x]
      
  #for i=0,n_elements(sion)-1 do begin
  #  m=where(uniq eq sion[i])
  #  a[m,m]=-1d0*sig[i]*b[i] + (-1d0*sig[i]*(1d0-b[i]))
  #endfor

  #for i=0,n_elements(dion)-1 do begin
  #  c=where(uniq eq dion[i])
  #  r=where(uniq eq sion[i])
  #  a[c,r]=sig[i]*b[i]
  #endfor
  
  return a
#####################################################################################  
def where(array,val):
  tog=0
  for x in range(len(array)):
    if array[x]==val:
      index=x
      tog=1
  if tog==0: 
    print "No match for ",val," in array "
    return
  return index
#####################################################################################     
def makearray(num,up,down):
  step=(up-down)/num
  array=[0]*int(num)
  for x in range(len(array)):
    val=step*float(x)
    array[x]=val+down
  array[0]=array[1]/float(2)
  return array
##################################################################################### 
def atoms2massfrac(n,solar,abu_array):
    mass_point=open('iso_u.dat','r')
    massf=mass_point.readlines()
    
    solar_n_point=open('lodders_2009_n_atoms.dat','r')
    solar_n=solar_n_point.readlines()
    
    solar_n=solar_n[1:len(solar_n)]
    
    for x in range(len(massf)):
      massf[x]=float(massf[x])
    
    nums=len(solar)
    index=[0]*len(solar)
    k=0
    for x in range(len(solar)):
      for y in range(len(abu_array)):
	if solar[x]==abu_array[y]: 
	  index[k]=y
	  k=k+1
    
    mass=[0]*len(solar)
    for x in range(len(solar)):
      mass[x]=massf[index[x]]
      
    summ=float(0)
    for x in range(len(solar_n)):
      solar_n[x]=float(solar_n[x])
      summ=summ+solar_n[x]*massf[x]
    
    for x in range(len(solar)):
      n[x]=n[x]*mass[x]/summ
      
    return n
#####################################################################################     
    
    
    

