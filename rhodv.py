from nedcalc1 import *
import os
from numpy import *
from calc_kcor import *
h=0.696
mgs=-19.39+5*log10(h)

f=open('impvalues','r')
lines=f.readlines()
f.close()
size=len(lines)
zimp=zeros(size-1)
ngal=zeros(size-1)
for i in range(1,size):
 line = lines[i].split('\n')[0]
 zimp[i-1] = float(line.split(',')[0])
 ngal[i-1] = float(line.split(',')[4])
 
f1=open('/home/adiman/Desktop/absorbers/folders','r')
folders=f1.readlines()
s1=len(folders)
f1.close()

for i in range(s1):
 system=folders[i].split('\n')[0]
 red = float(system.split('_')[1])
 sel = where(zimp==red)[0]
 if len(sel)!=0:
  if ngal[sel]>0:
   print('\n'+system)
   f=open('/home/adiman/Desktop/absorbers/'+system+'/galaxies/2mpc_1000kms_sorted.csv','r')
   lines1 = f.readlines()
   f.close()
   ng = 0
   f=open('/home/adiman/Desktop/absorbers/'+system+'/galaxies/rhodv','w')
   f.write('rho  rho/rvir dv  rho3d  rho3d/rvir')
   for j in range(len(lines1)):
    line = lines1[j].split('\n')[0].split()
    zgal = line[2]
    rho=float(line[-3])
    if rho<2000:
     if zgal!='*******':
      z = float(zgal)
      dv=abs(float(line[4]))
      if z>0:
       dl=cosmcalc(100*h,0.286,0.714,z)[1]
       color=float(line[-2])-float(line[-1])
       mpk=float(line[-2])+5-5*log10(1e+6*dl)
       k = calc_kcor('g',z,'g - r',color)
       mg = mpk-k
       lratio=10**(0.4*(mgs-mg))
       if lratio>0.2:
        rvirial=250*(lratio)**(0.2) #Prochaska et al.
        dl1=cosmcalc(100*h,0.286,0.714,red)[1]
        dlos=abs(dl1-dl)*1e+3
        rho3d = sqrt(dlos**2+rho**2)
        ratio=rho/rvirial
        ratio1=rho3d/rvirial
        f.write('\n'+str(float(round(rho,2)))+'  '+str(round(ratio,2))+'      '+str(float(round(dv,2)))+'  '+str(float(round(rho3d,2)))+'   '+str(round(ratio1,2)))    
   f.close()

f=open('impvalues1','w')
f.write('z,imp(kpc),ngals,rhobyrvir(min),|dv|')

sys = array([],dtype=str)
reds = array([])
for i in range(s1):
 system=folders[i].split('\n')[0]
 sys=append(sys,system)
 reds=append(reds,float(system.split('_')[1]))

print(reds)
for i in range(1,size):
 line = lines[i].split('\n')[0]
 z = float(line.split(',')[0])
 ngal = float(line.split(',')[4])
 print(z)
 if ngal==0:
  f.write('\n'+line)
 else:
  system=sys[reds==z]
  f1=open('/home/adiman/Desktop/absorbers/'+str(system[0])+'/galaxies/rhodv','r')
  lines1 = f1.readlines()
  ngal=len(lines1)-1
  if ngal!=0:
   rbyrvir=zeros(ngal)
   r=zeros(ngal)
   dvs=zeros(ngal)
   for j in range(1,len(lines1)):
    line = lines1[j].split('\n')[0].split()
    r[j-1]=float(line[0])
    rbyrvir[j-1]=float(line[1])
    dvs[j-1]=float(line[2])

   sel = where(rbyrvir==min(rbyrvir))[0]
   rbyrvirmin = min(rbyrvir) 
   dvs = dvs[sel] 
   dvmin = min(dvs)
   r = r[sel]
   rmin = float(r[dvs==dvmin])
   f.write('\n'+str(z)+','+str(rmin)+','+str(ngal)+','+str(rbyrvirmin)+','+str(dvmin))
  else:
   f.write('\n'+str(z)+',0,0,0,0') 
  
 
