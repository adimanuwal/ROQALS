#For determining line-of-sight velocity separation, k corrections, L/L* and virial radii of galaxies from SDSS spectroscopic database
from scipy.constants import c
from nedcalc1 import *#Ned Wright's Cosmology Calculator
from numpy import *
from calc_kcor import *#k corrections

sqlout = input('Enter the sql query output:')
red = float(input('Enter the reference redshift:'))#redshift about which the velocity separations are determined

h=0.696#Hubble parameter
mrs=-20.44+5*log10(h)#Mr* for z<0.16 based on Blanton et al. 2003

def relv(z1,z2):#relative los velocity separation in km/s
 dv=(c*1e-3)*((1+z2)**2-(1+z1)**2)/((1+z2)**2+(1+z1)**2)
 return(dv) 

#Sorting the galaxies acoording to their impact parameter
f=open(sqlout,'r')
lines=f.readlines()
f.close()

f=open('sorted.csv','w')#file with galaxies sorted by impact parameter
if len(lines)>2:#the galaxy data in a SDSS SQL query output begins after first two rows 
 ra=array([])#RA
 dec=array([])#Dec
 z=array([])#spectroscopic z
 zerr=array([])#error in spec z
 dv=array([])#los velocity separation
 dverr=array([])#error in dv
 angimp=array([])#angular separation in '
 linimp=array([])#linear separation in kpc
 g=array([])#g-band apparant magnitude
 r=array([])#r-band apparent magnitude
 for i in range(2,len(lines)):
   line=lines[i].split('\n')[0].split(',')
   ra=append(ra,round(float(line[1]),5))
   dec=append(dec,round(float(line[2]),5))
   z=append(z,round(float(line[5]),5))
   zerr=append(zerr,round(float(line[6]),5))
   dv=append(dv,round(relv(float(line[5]),red),4))
   num=(1+z[-1])**2-(1+red)**2
   den=(1+z[-1])**2+(1+red)**2
   numer=2*zerr[-1]
   dener=numer 
   dverr=append(dverr,round(abs(dv[-1]*sqrt((numer/num)**2+(dener/den)**2)),4))
   angimp=append(angimp,round(float(line[7]),4))
   imp=angimp[-1]*cosmcalc(100*h,0.286,0.714,red)[0]*60
   linimp=append(linimp,round(imp,4))
   g=append(g,round(float(line[10]),4))
   r=append(r,round(float(line[12]),4))

 #sorting
 sortind = argsort(linimp)
 ra=ra[sortind]
 dec=dec[sortind]
 z=z[sortind]
 zerr=zerr[sortind]
 dv=dv[sortind]
 dverr=dverr[sortind]
 angimp=angimp[sortind]
 linimp=linimp[sortind]
 g=g[sortind]
 r=r[sortind]
 for i in range(len(ra)):
   if i!=len(ra)-1:
    f.write(str(ra[i])+' '+str(dec[i])+' '+str(z[i])+' '+str(zerr[i])+' '+str(dv[i])+' '+str(dverr[i])+' '+str(angimp[i])+' '+str(linimp[i])+' '+str(g[i])+' '+str(r[i])+'\n')
   else: 
    f.write(str(ra[i])+' '+str(dec[i])+' '+str(z[i])+' '+str(zerr[i])+' '+str(dv[i])+' '+str(dverr[i])+' '+str(angimp[i])+' '+str(linimp[i])+' '+str(g[i])+' '+str(r[i]))
f.close()

lcut=0.01#Luminosity threshold for galaxies
#Saving the impact parameter, normalized impact parameter and velocity separation for galaxies above a luminosity cut
f=open('sorted.csv','r')
lines1 = f.readlines()
f.close()
f=open('rhodv','w')
f.write('rho  rho/rvir dv')
for j in range(len(lines1)):
   line = lines1[j].split('\n')[0].split()
   z = float(line[2])
   rho=float(line[-3])
   dv=abs(float(line[4]))
   if rho<1500 and dv < 1000 and z>0:#impact parameter < 1.5 Mpc and line-of-sight velocity separation < 1000 km/s
     dl=cosmcalc(100*h,0.286,0.714,z)[1]#luminosity distance
     color=float(line[-2])-float(line[-1])#(g-r) color
     mpk=float(line[-1])+5-5*log10(1e+6*dl)#Mr+k
     k = calc_kcor('r',z,'g - r',color)#k
     mr = mpk-k#k-corrected Mr
     lratio=10**(0.4*(mrs-mr))#L/L*
     if lratio>lcut:
       rvirial=250*(lratio)**(0.2) #R_vir (Prochaska et al. 2011)
       ratio=rho/rvirial
       f.write('\n'+str(float(round(rho,2)))+'  '+str(round(ratio,2))+'      '+str(float(round(dv,2))))    
f.close()
 
