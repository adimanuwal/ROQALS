#total pathlength and grating stats
import os
import numpy as np
import matplotlib.pyplot as plt
from COG import curve_of_growth
from scipy.constants import c
from scipy import stats
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
fh=0.416400
#lh=1215.6701
#f = open('cog_bpar13.out','r')
#lines = f.readlines()
#size = len(lines)
#logn = np.zeros(size)
#logw = np.zeros(size)

#for i in range(size):
 #line = lines[i].split('\n')[0].split()
 #logn[i] = np.log10(10**(floasot(line[0]))/(fh*lh))
 #logw[i] = np.log10(10**(float(line[1]))*lh)

#required constants
lr=1548.195
lr1=1550.770
#lr=lr1
fl=0.190800
fl1=0.095220
za=0.00167#for extragalactic source v>500 km/s

#COG
x,y = curve_of_growth(lr,fl, 13)
logn = np.log10(x/(fl*lr))
logw = np.log10(y*lr)

#Delta z
N=[14.9]#np.arange(12.8,14.9+0.1,0.1)
size=len(N)

#For npix corresponding to a given W
f1=open('civaodvalues','r')
lines1=f1.readlines()
s1=len(lines1)
cdet=s1-1

cz=np.zeros(cdet)
cdv=np.zeros(cdet)
ind=-1
for i in range(1,s1):
 line=lines1[i].split(',')
 ind+=1
 cz[ind]=float(line[1])
 cdv[ind]=float(line[5])-float(line[4])

f6=open('vpvalues','r')
lines6 = f6.readlines()
s6=len(lines6)
zvp=np.zeros(s6-1)
nciv=np.zeros(s6-1)
for i in range(1,s6):
 line=lines6[i].split(',')
 zvp[i-1]=float(line[0])
 nciv[i-1]=float(line[7])
czs = np.unique(cz)
size1 = len(czs)
cns = np.zeros(size1)
cdvs = np.zeros(size1)
for i in range(size1):
 cns[i] = round(np.log10(sum(10**nciv[zvp==czs[i]])),1)
 cdvs[i] = sum(cdv[cz==czs[i]])

tpl=open('/home/adiman/Desktop/Plots/corrvalues/dz_neg.txt','w')
nsigma=3 

#For mean width of all lines
os.chdir('/home/adiman/Desktop/absorbers/')
f1=open('folders','r')
folders=f1.readlines()
s1=len(folders)
#fetching the data
dva=np.array([])
for i in range(s1):
 system=folders[i].split('\n')[0]
 os.chdir('/home/adiman/Desktop/absorbers/'+system+'/vpfit')
 name=system.split('_')[0]
 z=system.split('_')[1]
 table=name+'_ltable.dat'
 f=open(table,'r')
 lines=f.readlines()
 s=len(lines)
 for j in range(1,s):
  line=lines[j].split()
  v1=float(line[12].split(',')[0].split('[')[1])
  v2=float(line[12].split(',')[1].split(']')[0])
  dva = np.append(dva,v2-v1)

#Contaminations from Danforth et al. 2016
cpath = '/home/adiman/danforth/'
f=open(cpath+'hlsp_igm_hst_cos_tab1.txt','r')
lines = f.readlines()
folders = []
names = []
for i in range(3,len(lines)):
 line = lines[i].split()
 folders.append(line[0])
 names.append(line[1])

#Sowgat's files
snames = []
sfiles = []
f=open('/home/adiman/Desktop/contfit/spectra/reqspec','r')
lines=f.readlines()

for i in range(len(lines)):
 line = lines[i].split('\n')[0].split()
 if len(line)==3:
  snames.append(line[0])
  sfiles.append(line[2])

print('mean dv:',np.mean(dva))
meandvp = 10
npixc=int(np.mean(dva)/meandvp)

os.chdir('/home/adiman/Desktop/contfit/')

f=open('/home/adiman/Desktop/Plots/corrvalues/modspec','r')
spectra=f.readlines()
s=len(spectra)

mdv = np.mean(cdvs)
for i in range(size):
 DZ=0
 #mdv = np.mean(cdvs)#np.mean(cdvs[np.where((cns>=N[i]-0.1) & (cns<=N[i]+0.1))[0]])
 print('\nAdding for log[N(CIV)] =',round(N[i],1))
 print('Mean dv (km/s):',mdv)
 nspectra = 0
 for l in range(s):
  name=spectra[l].split('\n')[0].split('.asc')[0]
  found=0
  dz=0
  for k in range(len(snames)):
   if snames[k]==name:
    fname = 'spectra/'+sfiles[k]
    found=1
  if found==0:
   fname = name+'_cf.dat'

  if os.path.exists('/home/adiman/Desktop/contfit/'+fname)==False:
   continue
  else:
   print(name)
  nspectra+=1
  zem=float(spectra[l].split('\n')[0].split('.asc')[1].split(',')[1])
  #Removing regions that should be avoided
  if 1<0:#len(folders[names==name])!=0:
   f = open(cpath+folders[names==name]+'/'+'hlsp_igm_hst_cos_'+folders[names==name]+'_g130m-g160m_v3_linelist.txt')
   linesc = f.readlines()
   wc = np.array([])
   ewc = np.array([])
   for j in range(74,len(linesc)):
    linec=linesc[j].split()
    if linec[2]!='1548"' or linec[1]=='NULL':
     wc = np.append(wc,float(linec[0]))
     ewc = np.append(ewc,float(linec[-11]))

   f1=open(fname,'r')
   lines1=f1.readlines()
   ws = np.array([])
   fluxs = np.array([])
   errs = np.array([])
   conts = np.array([])
   for j in range(len(lines1)):
    line1 = lines1[j].split()
    w=float(line1[0])
    dw = float(lines1[1].split()[0]) - float(lines1[0].split()[0])
    flag=0
    for k in range(len(wc)):
     if (w>wc[k]-0.5*npixc*dw) and (w<wc[k]+0.5*npixc*dw):
      flag=1
    if flag==0:
     ws = np.append(ws,w)
     fluxs = np.append(fluxs,float(line1[1]))
     errs = np.append(errs,float(line1[2]))
     conts = np.append(conts,float(line1[3]))
  else:
   f1=open(fname,'r')
   lines1=f1.readlines()
   ws = np.array([])
   fluxs = np.array([])
   errs = np.array([])
   conts = np.array([])
   for j in range(len(lines1)):
    line1 = lines1[j].split()
    w=float(line1[0])
    ws = np.append(ws,w)
    fluxs = np.append(fluxs,float(line1[1]))
    errs = np.append(errs,float(line1[2]))
    conts = np.append(conts,float(line1[3]))

  #plt.plot(ws,conts,'k-')
  #plt.plot(ws,errs,'r-')
  #plt.show()
  s1=len(ws)
  lmax=max(ws)
  zb=(lmax-lr1)/lr1#maximum z decided by 1550 coverage along with 1548

  j=0
  npix=0
  while j<s1-1: #in range(s1-1-5):
   w1=ws[j]
   z1=(w1-lr)/lr
   w2=ws[j+1]
   z2=(w2-lr)/lr
   vqso=(c*1e-3)*((1+z1)**2-(1+zem)**2)/((1+zem)**2+(1+z1)**2)#vqso>5000 km/s for non-association
   dvpix=(c*1e-3)*((1+z2)**2-(1+z1)**2)/((1+z2)**2+(1+z1)**2) 
   npix = mdv/dvpix
   dzi=0
   sbyn=0
   #W = 10**(np.interp(N[i],logn,logw))
   if w1>lr*(1+za) and w1<lr*(1+zb) and vqso<-5000 and np.round(w2-w1,2)<0.06 and j+int(npix/2.0)<=len(ws)-1:
    eqw = 0
    unc = 0
    flag = 0
    #unc=np.array([])
    #sel = np.array([],dtype=int)
    #for k in range(j-int(npix/2.0),j+int(npix/2.0)):
      #dl = ws[k+1]-ws[k]
      #sel = np.append(sel,k)
      #if np.round(dl,2)>=0.06:
       #flag=1
      #eqw=eqw+(1-fluxs[k]/conts[k])*dl
      #unc=np.sqrt(unc**2+(errs[k]/conts[k]*dl)**2)
     #sbyn = np.append(sbyn,(conts[k] - fluxs[k])/errs[k])
    #dl = w2-w1
    #eqw = (1-fluxs[j]/conts[j])*dl
    #unc = (errs[j]/conts[j]*dl)
    #unc = np.sqrt(sum(unc**2))
    sel = np.arange(j-int(npix/2.0),j+int(npix/2.0))
    sel1 = np.arange(j-int(npix/2.0),j+int(npix/2.0)+1)
    wsel = ws[sel]
    dls = np.diff(ws[sel1])
    eqws = (1-fluxs[sel]/conts[sel])*dls
    uncs = errs[sel]/conts[sel]*dls
    eqw = sum(eqws)
    unc = np.sqrt(sum(uncs**2))
    sbyn = eqw/unc
    #print(unc)
    #if unc!=0:
     #sbyn = np.append(sbyn,eqw/unc)
    #if sbyn<0:
     #plt.step(ws[sel],fluxs[sel]/conts[k],'b-')
     #plt.ylabel('Normalized Flux',fontsize=15)
     #plt.xlabel(r'$\lambda$',fontsize=15)
     #plt.show()
    #if abs(sbyn)>0 and unc!=0:# and flag==0:
    if unc!=0:# and flag==0:
     #sbyn = eqw/unc 
     #Wlim = nsigma*np.sqrt(np.ceil(npix))*(w2-w1)/(sbyn*(1+z1))#Hellsten+98
     Wlim = 3*unc/(1+z1)#Burchett+15
     #Wlim = 3*dl/sbyn#Danforth+08
     #lgN=np.log10(Wlim)-np.log10(lr1**2*fl1)+20.053#Jill Bechtold
     lgN = np.interp(np.log10(Wlim),logw,logn)
     if lgN<=N[i]:
      dzi=(z2-z1)
      dz+=dzi
   j+=1
  DZ+=dz
  print('Current DZ:',DZ)
 print('No. of spectra:',nspectra)
 print('DZ=',round(DZ,2))
 tpl.write(str(round(N[i],1))+'  '+str(DZ)+'\n')

 


   
      

    
