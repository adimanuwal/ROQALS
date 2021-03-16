#Absorption Line Analysis based on the Apparent Optical Depth (AOD) method
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
import os
cval=c
#function for converting to scientific notation
def sci(n):
 a='%e' % n
 sa=str(a)
 spa=sa.split('e')
 m=round(float(spa[0]),2)
 if np.shape(spa)[0]>1:
  scn=str(m)+'e'+spa[1]
 else:
  scn=str(m)+'e0'
 return scn
#function for step filling
def fill_between_steps(x, y1, y2, h_align='right', ax=None, **kwargs):
    ''' Fills a hole in matplotlib: Fill_between for step plots.

    Parameters :
    ------------

    x : array-like
        Array/vector of index values. These are assumed to be equally-spaced.
        If not, the result will probably look weird...
    y1 : array-like
        Array/vector of values to be filled under.
    y2 : array-Like
        Array/vector or bottom values for filled area. Default is 0.

    **kwargs will be passed to the matplotlib fill_between() function.

    '''
    # If no Axes opject given, grab the current one:
    if ax is None:
        ax = plt.gca()
    # First, duplicate the x values
    xx = x.repeat(2)[1:]
    # Now: the average x binwidth
    xstep = (x[1:] - x[:-1]).mean()
    # Now: add one step at end of row.
    xx = np.append(xx, xx.max() + xstep)

    # Make it possible to change step alignment.
    if h_align == 'mid':
        xx -= xstep / 2.
    elif h_align == 'right':
        xx -= xstep

    # Also, duplicate each y coordinate in both arrays
    y1 = y1.repeat(2)
    if type(y2) == np.ndarray:
        y2 = y2.repeat(2)

    # now to the plotting part:
    return ax.fill_between(xx, y1, y2=y2, **kwargs)

#For vertical lines
def drawline(x):
    a = np.ones(10)*x
    b = np.linspace(0,1.5,10)
    plt.plot(a,b)

#For equivalent width calculation
def wandn(v1,v2):
     print('v1:',int(round(v1)),'km/s')
     print('v2:',int(round(v2)),'km/s')
     no=0
     for i in range(s):
      if v[i]>=v1 and v[i]<=v2:
       no=no+1
     m=np.zeros(no)
     n=np.zeros(no)
     ind=0
     eqw=0
     unc=0
     for j in range(s):
      if v[j]>=v1 and v[j]<=v2:
       m[ind]=v[j]
       n[ind]=nflux[j]
       ind+=1
       dl=(wave[j+1]-wave[j])*1e+3
       eqw=eqw+(1-nflux[j])*dl
       unc=np.sqrt(unc**2+(nerr[j]*dl)**2)
     eqwr=int(eqw/(1+z2))
     uncr=int(round(abs(unc/(1+z2))))
     flag=0
     if eqwr<3*uncr:
      flag=1
      print('Non-detection!')
      print('Determining 3 sigma upper limit of W ...')
      uncr1=uncr
      eqw=3*unc
      eqwr=3*uncr
      unc=0
      uncr=0
     print('Observed Equivalent Width (mAngs):',int(eqw),' +- ',int(round(unc)))
     print('Rest Equivalent Width (mAngs):',eqwr,' +- ',uncr)
     datf.write(str(eqwr)+' +- '+str(uncr)+'    ')
     #N calculation
     ind=0
     ln=0
     ld=0
     punc=0
     if flag==1:
      print('Determining upper limit of N using the COG ...') 
      acd=3*uncr1/((twav/(1+z2))**2*fl)*10**(20.053-3) #by Jill Bechtold 
      lunc=0
     else:
      acd=0
      for j in range(s):
       if v[j]>=v1 and v[j]<=v2:
        m[ind]=v[j]
        n[ind]=nflux[j]
        ind+=1
        tau=-np.log(nflux[j])
        ad=3.768e+14*(fl*twav/(1+z2))**(-1)*(tau)
        dv=v[j+1]-v[j]
        acd=acd+ad*dv
        punc=np.sqrt(punc**2+(3.768e+14*(fl*twav/(1+z2))**(-1)*(-nerr[j]*dv/nflux[j]))**2)
      lunc=punc/(np.log(10)*acd)#base10
     print('Apparant Column Density (N(/sq.cm)):',sci(acd),' +- ',sci(abs(punc)))
     print('Apparant Column Density (log(N)(/sq.cm)):',str(round(np.log10(abs(acd)),2)),' +- ',str(round(abs(lunc),2)))
     datf.write(str(sci(acd))+' +- '+str(sci(abs(punc)))+'    '+str(round(np.log10(abs(acd)),2))+' +- '+str(round(abs(lunc),2))+'     ')
     fill_between_steps(m,n,1,color='red')
#velocity width
def kmatics(v1,v2):
     #print('v1:',v1,'km/s'
     #print('v2:',v2,'km/s'
     intau=0
     meanv=0
     omvs=0
     interr=0
     mverr=0
     omverr=0
     v1n=0
     v2n=0
     #for integrated OD
     for j in range(s):
      if v[j]>=v1 and v[j]<=v2:
       tau=-np.log(nflux[j])
       dv=v[j+1]-v[j]
       intau=intau+tau*dv
       interr=np.sqrt(interr**2+(nerr[j]*dv/nflux[j])**2)
     #for mean
     merr=0
     for j in range(s):
       if v[j]>=v1 and v[j]<=v2:
        tau=-np.log(nflux[j])
        dv=v[j+1]-v[j]
        meanv=meanv+tau*v[j]*dv
        merr=np.sqrt(merr**2+(nerr[j]*v[j]*dv/nflux[j])**2)
     mverr=abs((merr/intau)-interr*meanv/intau**2)
     meanv=meanv/intau
     #for omverr
     omv1err=0
     for j in range(s):
       if v[j]>=v1 and v[j]<=v2:
        tau=-np.log(nflux[j])
        dv=v[j+1]-v[j]
        omvs=omvs+(tau*((v[j]-meanv)**2)*dv)
        omv1err=np.sqrt(omv1err**2+((2*(v[j]-meanv)*mverr*tau+(v[j]-meanv)**2*nerr[j]/nflux[j])*dv)**2)
     if omvs<0:
      print('Negative value of omvs!') 
     omv=np.sqrt(omvs/intau)
     omverr=abs((omv1err/intau-interr*omvs/intau**2)/(2*omv))
     #for dv90
     t=0
     for j in range(s):
       if v[j]>=v1 and v[j]<=v2 and t<=(.05*intau):
        tau=-np.log(nflux[j])
        dv=v[j+1]-v[j]
        t=t+tau*dv
        v1n=v[j]
     t=0
     for j in range(s):
       if v[j]>=v1 and v[j]<=v2 and t<=(.95*intau):
        tau=-np.log(nflux[j])
        dv=v[j+1]-v[j]
        t=t+tau*dv
        v2n=v[j]
    
     dv90=v2n-v1n
     print('Integrated tau:',round(intau,2),'+-',round(interr,2))
     print('Mean v(km/s):',round(meanv,2),'+-',round(mverr,2))
     print('Omega v(km/s):',round(omv,2),'+-',round(omverr,2))
     print('dv90(km/s):',round(dv90,2))
     velar=np.zeros(2)
     velar[0]=int(round(v1))
     velar[1]=int(round(v2))
     datf.write('['+str(velar[0])+','+str(velar[1])+']'+'   ')
     datf.write(str(round(intau,2))+' +- '+str(round(interr,2))+'  ')
     datf.write(str(round(meanv,2))+' +- '+str(round(mverr,2))+'   ')
     datf.write(str(round(omv,2))+' +- '+str(round(omverr,2))+'   ')
     datf.write(str(round(dv90,2))+'\n')
#For interfacing
def onclick(event):
    thisline = event.artist
    xdata = thisline.get_xdata()
    ydata = thisline.get_ydata()
    ind = event.ind
    x=float(xdata[ind])
    y=float(ydata[ind])
    global v1
    global v2
    global c
    if c==0:
     v1=x
     c=1
    if x!=v1 and x>v1:
     v2=x
    drawline(x)
    #For coloring the region
    if v1<v2 and v2!=-1:
     wandn(v1,v2)
     kmatics(v1,v2)
    fig.canvas.draw()

path='/home/adiman/projectdata/COS.Legacy.Data/BIN.DATA/MODIFIED/contfit/'#path to the continuum fitted files
os.chdir(path)
#atlist=input('Enter the ion list file:')
atlist='modatoms.dat'#files with the atomic transition data
f=open(atlist,'r')
rows=f.readlines()
n=len(rows)
new='y'
while new=='y':
 fname=input('Spectrum(continuum fitted):')
 fsp=fname.split('_')
 datf=open(fsp[0]+'_ltable.dat','w')
 z2=float(input('z(absorber):'))
 print('Observe the system plot.')
 datf.write('Line\t   W(mAngs)    N(/sq.cm)\t       log(N)(/sq.cm)    [-v, +v]       int.tau\t       meanv(km/s)     omv(km/s)       dv90(km/s)\n')
 f=open(fname)
 lines=f.readlines()
 s=len(lines)
 wave=np.zeros(s)
 flux=np.zeros(s)
 err=np.zeros(s)
 cont=np.zeros(s)
 nflux=np.zeros(s)
 nerr=np.zeros(s)
 v=np.zeros(s)
 ind=-1
 #data extraction and normalized value calculation
 for line in lines:
   ind=ind+1
   el=line.split()
   wave[ind]=float(el[0])
   flux[ind]=float(el[1])
   err[ind]=float(el[2])
   cont[ind]=float(el[3])
   nflux[ind]=flux[ind]/cont[ind]
   nerr[ind]=err[ind]/cont[ind]
   if nflux[ind]<nerr[ind] and nflux[ind]<=0:
    nflux[ind]=nerr[ind]
 ch='y'
 while ch=='y':
  fig = plt.figure(figsize=(10,5))
  pos  = [0.08, 0.13, 0.9, 0.85] ; ax  = fig.add_axes(pos)
  trans=input('Transition:')
  datf.write(trans+'   ')
  const=trans.split()
  e0=const[0]
  e1=const[1]
  e2=const[2]
  check=0
  for i in range(n):
   r=rows[i].split()
   r2=r[3].split('.')
   if e0==r[0] and e1==r[1] and e2==r2[0]:
    print('Transition found in the list.')
    twav=float(r[3])*(1+z2)
    fl=float(r[4])
    check=1
  if check==0:
   print('Transition is not in the list!')
  if twav<min(wave) or twav>max(wave):
   print('The required observed wavelength is not covered in the given file !')
  if twav>=min(wave) and twav<=max(wave) and check==1:
      linw=twav/(1+z2)
      print('Rest transition wavelength(Angs):',linw)
      for m in range(s):
       z1=(wave[m]-linw)/linw
       v[m]=(cval*1e-3)*((1+z1)**2-(1+z2)**2)/((1+z1)**2+(1+z2)**2)
      ax.text(-235,0.2,trans,color='red',size=20)
      ax.text(150,0.2,round(linw*(1+z2),4),color='red',size=20)
      line, = ax.step(v,nflux,'black',picker=5)
      ax.step(v,nerr,'blue')
      #plt.tick_params(axis='y',which='major',pad=15)
      #plotting faint continuum and vertical
      ncont=np.ones(s)
      plt.plot(v,ncont,'--',color='black')
      verty=np.linspace(0,1.5,10)
      vertx=np.zeros(10)
      plt.plot(vertx,verty,'-.',color='black')
      plt.xlim(-300,300)
      plt.ylim(0,1.5)
      plt.xlabel('Velocity (km/s)',size=18)
      plt.ylabel('Normalized Flux',size=18)
      plt.xticks(size=15)
      plt.yticks(size=15)
      #line, = ax.plot(np.random.rand(100), 'o', picker=5)  # 5 points tolerance
      v1=-1
      v2=-1
      click=input('Velocity limits by clicking on the figure?(y/n):')
      if click=='y':
       c=0#flag
       print('Plotting',trans,'..')
       fig.canvas.mpl_connect('pick_event', onclick)
       plt.show()
      if click=='n':
       ax.set_title('Selected Region', size=20)
       v1=float(input('v1(km/s):'))
       drawline(v1)
       v2=float(input('v2(km/s):'))
       drawline(v2)
       wandn(v1,v2)
       kmatics(v1,v2)
       plt.show()
  ch=input('Another transition in the same file?(y/n):')
 print('Output file: '+fsp[0]+'_ltable.dat')
 new=input('New file?(y/n):')
 
