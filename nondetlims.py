from nedcalc1 import *
import os
from numpy import *
from calc_kcor import *
from scipy.constants import c

h=0.696
lr=1548.195
fl=0.190800
mgs=-19.39+5*log10(h)

def relv(z1,z2):
 dv=(c*1e-3)*((1+z2)**2-(1+z1)**2)/((1+z2)**2+(1+z1)**2)
 return(dv) 

for i in range(40):
 if os.path.exists('/home/adiman/Desktop/Plots/corrvalues/nondetgals/2mpc_1000kms_sorted_'+str(i+1)+'.csv')==True:
   f=open('/home/adiman/Desktop/Plots/corrvalues/nondetgals/2mpc_1000kms_sorted_'+str(i+1)+'.csv','r')
   lines1=f.readlines()
   f.close()
   f=open('/home/adiman/Desktop/Plots/corrvalues/nondetgals/rhodv_'+str(i+1),'w')
   print('rhodv_'+str(i+1))
   f.write('rho  rho/rvir dv z')
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
        ratio=rho/rvirial
        f.write('\n'+str(float(round(rho,2)))+'  '+str(round(ratio,2))+'      '+str(round(dv,2))+' '+str(z))    
   f.close()


f1 = open('nondetsights','r')
lines = f1.readlines()
f1.close()

#Sowgat's files
snames = []
sfiles = []
f=open('/home/adiman/Desktop/contfit/spectra/reqspec','r')
slines=f.readlines()

for i in range(len(slines)):
 line = slines[i].split('\n')[0].split()
 if len(line)==3:
  snames.append(line[0])
  sfiles.append(line[2])

f=open('nondetimpvalues','w')
f.write('QSO,z,imp(kpc),ngals,rhobyrvir(min),W,Werr,lgN,lgNerr')

for i in range(40):
 if os.path.exists('/home/adiman/Desktop/Plots/corrvalues/nondetgals/rhodv_'+str(i+1))==True:
  f1=open('/home/adiman/Desktop/Plots/corrvalues/nondetgals/rhodv_'+str(i+1),'r')
  lines1 = f1.readlines()
  ngal=len(lines1)-1
  if ngal!=0:
   for j in range(1,len(lines1)):
    line = lines1[j].split('\n')[0].split()
    r=float(line[0])
    rbyrvir=float(line[1])
    dv=float(line[2])
    z=float(line[3])

    if rbyrvir<=1.4:#if it is close by
     #W and N limits
     name = lines[i].split('\n')[0].split()[0]
     found=0
     for k in range(len(snames)):
      if snames[k]==name:
       fname = '/home/adiman/Desktop/contfit/spectra/'+sfiles[k]
       found=1
     if found==0:
      fname = '/home/adiman/Desktop/contfit/'+name+'_cf.dat'
     fc=open(fname,'r')
     linesc=fc.readlines()
     ws = array([])
     fluxs = array([])
     errs = array([])
     conts = array([])
     for k in range(len(linesc)):
      line1 = linesc[k].split()
      w=float(line1[0])
      ws = append(ws,w)
      fluxs = append(fluxs,float(line1[1]))
      errs = append(errs,float(line1[2]))
      conts = append(conts,float(line1[3]))
     zs = (ws-lr)/lr
     vs = relv(z,zs)
     sel = where((vs>-105/2.0) & (vs<105/2.0))[0]#select 105 km/s around the galaxy redshift
     vs = vs[append(sel,max(sel)+1)]
     dls = diff(ws[append(sel,max(sel)+1)])*1e+3
     fluxs = fluxs[sel]
     errs = errs[sel]
     conts = conts[sel]
     nflux = fluxs/conts
     nerr = errs/conts

     eqws = (1-nflux)*dls
     uncs = nerr*dls
     eqw = sum(eqws)
     unc = sqrt(sum(uncs**2))
     eqwr=int(eqw/(1+z))
     uncr=int(round(abs(unc/(1+z))))
     if eqwr<3*uncr:
      acd=round(log10(3*uncr/(lr**2*fl)*10**(20.053-3)),2) #by Jill Bechtold 
      lunc=0
      eqwr=int(3*uncr)
      uncr=0
      print('ND')
      print(acd,lunc)
     else:
      acd=0
      punc=0
      for k in range(len(nflux)):
       if nflux[k]>0:
        tau=-log(nflux[k])
        ad=3.768e+14*(fl*lr)**(-1)*(tau)
        dv=vs[k+1]-vs[k]
        acd=acd+ad*dv
        punc=sqrt(punc**2+(3.768e+14*(fl*lr)**(-1)*(-nerr[k]*dv/nflux[k]))**2)  
      lunc=round(punc/(log(10)*acd),2)#base10    
      acd = round(log10(acd),2)
      print('D')
      print(acd,lunc) 
     f.write('\n'+name+','+str(float(z))+','+str(float(r))+','+str(ngal)+','+str(float(rbyrvir))+','+str(int(eqwr))+','+str(int(uncr))+','+str(acd)+','+str(lunc))
   #else:
    #f.write('\n'+str(z)+',0,0,0,0') 
  
 
