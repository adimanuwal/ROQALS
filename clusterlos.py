#Finding COS sightlines that pass within 2 Mpc in projected separation wrt to BCGs in Abell Clusters 
from nedcalc1 import *
from numpy import *
f0 = open('modspec','r')
lines0 = f0.readlines()
f = open('abell_new','r')
lines = f.readlines()
f1 = open('qsoradec','r')
lines1 = f1.readlines()
f2 = open('clustercover','w')
cf = pi/180.0
f2.write('BCG               z       N(<1Mpc) N_sight  Distance [Mpc]  Sightline\n')
n1 = 0
n2 = 0
n3 = 0
for i in range(28,len(lines)):
  print('Cluster '+str(i-27)+' ...')
  line = lines[i].split()
  z = float(line[4])
  zsp = float(line[5])
  nmem = int(line[7])
  zsel = zsp
  if zsel ==-1:
   zsel = z - 0.03
  ra1 = float(line[2])*cf
  dec1 = float(line[3])*cf
  sf = cosmcalc(69.6,0.286,0.714,zsel)#scale factor in kpc/"
  count=0
  sname = array([],dtype=str)
  ds=array([])
  for j in range(len(lines1)):
   line1 = lines1[j].split()
   line0 = lines0[j].split('\n')[0]
   zqso = float(line0.split(',')[1])
   if zqso>zsel:
    ra2 = float(line1[1])*cf
    dec2 = float(line1[2])*cf
    theta = 180/pi*arccos(sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*cos(ra1-ra2))*3600
    d = sf*theta*1e-3#projected distance in Mpc
    if d<=2:#if the sightline is within 2 Mpc of the BCG
     count+=1
     sname = append(sname,line1[0])
     print('yes!: '+sname[-1]+','+str(round(d,1))+' Mpc')
    #if d<dmin:
     ds=append(ds,d)
    if d<1:  
     n1+=1
    if d>1 and d<1.5:
     n2+=1
    if d>1.5 and d<2:
     n3+=1 
  f2.write('\n'+line[1]+'  '+str(zsel)+'   '+str(nmem)+'       '+str(count)+'        '+str(ds)+'             '+str(sname))
f2.close()
   
print('Sightlines at d <1 Mpc:',n1)
print('Sightlines at 1<d<1.5 Mpc:',n2)
print('Sightlines at 1.5<d1<2 Mpc:',n3)
  
  
