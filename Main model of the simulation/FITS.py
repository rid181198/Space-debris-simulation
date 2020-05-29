# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 01:00:21 2020

@author: Admin
"""

from astropy.io import fits
from astropy import wcs
import numpy as np
import random
import Noises as ns
import final2 as f
import config as c
import matplotlib.pyplot as plt 
import math as m 
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
num=105
filename = 'C:/Users/Admin/Desktop/'+ str(num)+'.fits'
# making an image with 1024 X 1024 pixels with 16 bit

# data = np.zeros((1024, 1024), dtype=np.uint16)
# hdu = fits.PrimaryHDU(data=data)
# hdu.writeto(filename)
# hdul = fits.open(filename)

# opening the file as list
with fits.open(filename) as hdul:
    data = hdul[0].data


    
# #i = rows, altitude
# #j = column, azimuth



#field of view second method
#fn is f# number
fnm = 13
#magnification
M=1
fn=fnm*M
#diameter of telescope in mm
D = 1200
no_pixel = 1024
#pixel size in micron
spixel = 13
#plate scale in arc sec/mm
plate_scale = 206265/(D*fn)
#per pixel
pp = plate_scale*spixel/1000
#field of view
fov2 = no_pixel*pp*0.000277778
print(fov2)
print(fov2/0.000277778)



# #size of ccd chip
# #vertical size 30 mm
# ccdv = 13.5
# #horizontal size 30 mm
# ccdh = 30
# #focal length of telescope mm
# focal = 6000
# #field of view
# fov = 2*m.atan(ccdv/(2*focal))
# #size of field in km
# d=420
# fov2 = d*fov
# fov2 = fov*180/np.pi
# print(fov2)
# #fov2/pixel



pixel = 1024
P = fov2/pixel



#size of the satellite in km
size = 0.025
#distance in km
d =420
#angle
alpha = size*180/(d*np.pi)

pcoverage = alpha/P
# print(pcoverage)
s=4


#i = columns
#j = rows


f.alt_az.position(c.data[0],c.data[1],c.data[2],c.data[3],c.data[4],c.data[5]+4,c.data[6],c.data[7])
# distance = np.linalg.norm(sat_obs_vector)
alt0 = f.alt_az.elevation
az0 = f.alt_az.azimuth
# print(alt0)
# print(az0)
#field of view
alt1 = alt0 + fov2/2
alt2 = alt0 - fov2/2
az1 = az0 + fov2/2
az2 = az0 - fov2/2

#pos
posalt=[]
posaz=[]
ALT=[]
RA=[]
DEC=[]
# print(alt1,alt2)
# print(az1,az2)

#list of position and time
tl=[]
# maxp = int(np.amax(photon_l))

adc_l=[]
for t in range(0,1000,1):
    t=t/100
    f.alt_az.position(c.data[0],c.data[1],c.data[2],c.data[3],c.data[4],c.data[5]+t,c.data[6],c.data[7])
    f.ra_dec.position(c.data[0],c.data[1],c.data[2],c.data[3],c.data[4],c.data[5]+t)
    
    # distance = np.linalg.norm(sat_obs_vector)
    if (alt2<= f.alt_az.elevation <=alt1 and az2<= f.alt_az.azimuth <=az1):
        tl.append(t)
        posalt.append(f.alt_az.elevation)
        posaz.append(f.alt_az.azimuth)
        
        ALT.append(f.alt_az.altitude)
        RA.append(f.ra_dec.RA)
        DEC.append(f.ra_dec.declination)
        f.photons.get(t,100,0.8,0.0070)
        adc = (f.photons.num)*65535/(100000*3.14*5*5)
        
         
        adc_l.append(adc)
    # posalt.append(f.alt_az.elevation)
    # posaz.append(f.alt_az.azimuth)
#time of flight in the field of view

dt = tl[len(tl)-1] - tl[0]
print(dt)
#trails
alto=[]


for j in range(0,pixel,1):
    for alt in posalt:
        d1 = (alt1-(j*fov2)/pixel) -alt
        if (abs(d1) <= 0.001 and j not in alto):
            # print(abs(d1))
            alto.append(j)
           
azo=[]
for i in range(0,pixel,1):
    for az in posaz:
        d2 = (az2 + (i*fov2)/pixel) -az
        
        if (abs(d2) <= 0.001 and i not in azo):
            # print(abs(d2))
            azo.append(i)
            
          
nazo=[]
nalto=[]


if (len(azo)>=len(alto)):
    l = len(azo)/len(alto)
    nalto = alto[:]
    for j in np.arange(0,float(len(azo)),l):
        nazo.append(azo[int(j)])
    dp=nalto[2]-nalto[0]
if (len(azo)<len(alto)):
    l = len(alto)/len(azo)
    nazo = azo[:]
    for j in np.arange(0,float(len(alto)),l):
        nalto.append(alto[int(j)])
    dp = nazo[1]-nazo[0] 
    

  




slopel =[]
cl=[]
for i in range(0,len(nazo)-1):
    slope = (nalto[i+1]-nalto[i])/(nazo[i+1]-nazo[i])
    slopel.append(slope)
    
    C = nalto[i] - (slope*nazo[i])
    cl.append(C)


nazo2l=[]
nalto2l=[]
length = len(nazo)
for i in range(0,length-1):
    for x in range(nazo[i],nazo[i+1],1):
        for y in range(nalto[i],nalto[i+1],1):
            
            if (y == int(slopel[i]*x + cl[i])):
                nazo2l.append(x)
                nalto2l.append(y)



#for per pixel satllite time
range_azo = max(nazo2l) - min(nazo2l)
range_alto = max(nalto2l) - min(nalto2l)
cov = np.power(range_azo*range_azo + range_alto*range_alto,1/2)

pdt = dt/cov
print(pdt)
def tempcov(x,y):
    rangex = x - min(nazo2l)
    rangey = y - min(nalto2l)
    cov = np.power(rangex*rangex + rangey*rangey,1/2)
    return cov



    
#satellite imaging

adu =0
#argument as follows pixels, exposure time, quantum efficiency, area of pixel in micron, photon flux 
ns.noise.get(1024,dt,60,169,1)
datasky = np.random.poisson(ns.noise.sk,1024*1024)   

datadark = np.random.poisson(ns.noise.nd,1024*1024)


fwhm = 1.028*s/1.22
# fwhm = 1.028*s/1.22
sigma = fwhm/2.355
for i in range(0,1024,1):
        for j in range(0,1024,1): 
            #noises : dark noise, sky background, readout noise,                              
            data[i,j]=(datadark[j*i]*65535/(100000)) + (datasky[j*i]*65535/(100000))+ abs((np.random.normal(0, ns.noise.nr)*65535/100000))  + abs((np.random.normal(0, ns.noise.nw)*65535/100000))    +abs((np.random.normal(0, ns.noise.nf)*65535/100000))    



count = 0
for (z,x) in zip(nazo2l,nalto2l):                

    s=5
    #center of the disc
    r = x
    l =z
    
    
    K=s*s
  
    count = count +1
    f.photons.get(tempcov(l,r)*pdt,35,0.8,pdt*12)
   
    # ns.noise.get(1024,dt,60,169,1)
    adc = np.random.poisson(f.photons.num*65535/(100000*3.14*5*5),1500)
  
    # adc = (adc)*65535/(100000*3.14*5*5) 
    for i in range(r-2*s,r+2*s,1):
        for j in range(l-2*s,l+2*s,1):
              if (0<=r-2*s <= pixel and 0<=r+2*s <= pixel and 0<=l-2*s <= pixel and 0<=l+2*s <= pixel):
                 
                  if ( 0 <=np.power(i-r,2) + np.power(j-l,2) <= K+K):
                     
                      adu =  data[i,j] + int(((adc[750])) * (m.exp(-1*( (np.power(i-r,2) + np.power(j-l,2))/(2*np.power(sigma,2))) )))
                      # print(((adc[i])*65535/(100000*3.14*5*5)))
                      if (adu > max(adc)):
                          adu= max(adc)
                      data[i,j]=adu
      
            
                     
                        
             
        
#making list of images or objects
hdu1=fits.PrimaryHDU(data=data)        
# hdu2 = fits.ImageHDU(data=data)
new_hdul = fits.HDUList([hdu1])

#writing details
hdr = new_hdul[0].header
hdr['TYPE']=('Satellite Trails' , 'The purpose')

hdr['DATE'] = (str(c.data[0])+'/'+str(c.data[1]) +'/'+str(c.data[2]) +' '+ str(c.data[3]) + ' hh : ' + str(c.data[4]) + ' mm : ' + str(c.data[5]+4) + ' ss'   ,'Time and Date in universal time')    
hdr['FOV'] = (str(round(fov2,2)) + ' X '+ str(round(fov2,2))+' degrees', 'Field of view of CCD chip')
hdr['ALTITUDe'] = (str(round(alt2,2))+' X '+str(round(alt1,2))+' degrees', 'Altitude range of the field of view')
hdr['AZIMUTH'] = (str(round(az2,2))+' X '+str(round(az1,2))+' degrees','Azimuth range of field of view')
hdr['DISTANCE'] = (str(d)+' km','Distance between observer and satellite')
hdr['TIME'] = (str(round(dt,2)) + ' sec','Time of passing through field of view')

hdu1 = fits.PrimaryHDU(data=data,header=hdr)
# hdu2 = fits.ImageHDU(data=data,header=hdr)
new_hdul = fits.HDUList([hdu1])
#overwrite 
new_hdul.writeto(filename,overwrite=True)



image_data = fits.getdata(filename)
image_header = fits.getheader(filename,0)
print(image_header)


plt.imshow(image_data, cmap='gray',vmin = 0, vmax = 65535)
plt.colorbar()


             
# hdu1=fits.PrimaryHDU(data=data)        

# new_hdul = fits.HDUList([hdu1]) 
# new_hdul.writeto(filename,overwrite=True)
# image_data = fits.getdata(filename)



plt.imshow(image_data, cmap='gray',vmin = 0, vmax = 65535)
plt.colorbar()













#inclination
v=(90+abs(DEC[4]))*np.pi/180
u=RA[4]*np.pi/180
x=np.cos(u)*np.sin(v)
y=np.sin(u)*np.sin(v)
z=np.cos(v)

v=(90+abs(DEC[10]))*np.pi/180
u=RA[10]*np.pi/180
X=np.cos(u)*np.sin(v)
Y=np.sin(u)*np.sin(v)
Z=np.cos(v)

vec1 = (X-x,Y-y,Z-z)
vec1 = vec1/(np.linalg.norm(vec1))

v=90*np.pi/180
u=RA[4]*np.pi/180
x=np.cos(u)*np.sin(v)
y=np.sin(u)*np.sin(v)
z=np.cos(v)

v=90*np.pi/180
u=RA[10]*np.pi/180
X=np.cos(u)*np.sin(v)
Y=np.sin(u)*np.sin(v)
Z=np.cos(v)

vec2 = (X-x,Y-y,Z-z)
vec2 = vec2/(np.linalg.norm(vec2))

print('Inclination : ' + str(m.acos(np.dot(vec1,vec2))*180/np.pi) + ' degree')














#velocity and distance
v=(90+abs(DEC[0]))*np.pi/180
u=RA[0]*np.pi/180
x=np.cos(u)*np.sin(v)
y=np.sin(u)*np.sin(v)
z=np.cos(v)

v=(90+abs(DEC[-1]))*np.pi/180
u=RA[-1]*np.pi/180
X=np.cos(u)*np.sin(v)
Y=np.sin(u)*np.sin(v)
Z=np.cos(v)

vec1 = (X-x,Y-y,Z-z)
# cubic root of mass of the earth and gravitational constant
vel = np.power((np.linalg.norm(vec1)/dt)*6.674*np.power(10,-20.0)*5.972*np.power(10,24.0),1/3)
print('Velocity : ' + str(vel) + ' km/s')

dist = vel/(np.linalg.norm(vec1)/dt)
print('Distance : ' + str(dist ) + ' km')













#areal velocity of satellite and areal velocity
v=(90+abs(DEC[0]))*np.pi/180
u=RA[0]*np.pi/180
x=np.cos(u)*np.sin(v)
y=np.sin(u)*np.sin(v)
z=np.cos(v)


v=(90+abs(DEC[-1]))*np.pi/180
u=RA[-1]*np.pi/180
X=np.cos(u)*np.sin(v)
Y=np.sin(u)*np.sin(v)
Z=np.cos(v)

vec1 = (X-x,Y-y,Z-z)
# cubic root of mass of the earth and gravitational constant
vel = np.power((np.linalg.norm(vec1)/dt)*6.674*np.power(10,-20.0)*5.972*np.power(10,24.0),1/3)
length = vel*dt


arealv = (dist)*length/(2*dt)









i=3

j=i+1
#eccentricity
v=(90+abs(DEC[i]))*np.pi/180
u=RA[i]*np.pi/180
x=np.cos(u)*np.sin(v)
y=np.sin(u)*np.sin(v)
z=np.cos(v)


v=(90+abs(DEC[j]))*np.pi/180
u=RA[j]*np.pi/180
X=np.cos(u)*np.sin(v)
Y=np.sin(u)*np.sin(v)
Z=np.cos(v)



velv = (x-X,y-Y,z-Z)                 
posv = (x,y,z)

velvm = np.linalg.norm(velv)
posvm = np.linalg.norm(posv)

vel = velv/velvm
pos = posv/posvm

angle = m.acos(np.dot(vel,pos))

print(angle*180/np.pi)

vec1 = (X-x,Y-y,Z-z)
dt = tl[j] - tl[i]
# cubic root of mass of the earth and gravitational constant
velocity = np.power((np.linalg.norm(vec1)/dt)*6.674*np.power(10,-20.0)*5.972*np.power(10,24.0),1/3)
dist = velocity/(np.linalg.norm(vec1)/dt)


G = 6.674*np.power(10,-11.0)
M = 5.972*np.power(10,24.0)

E = np.power(np.power((dist*velocity*velocity/(G*M)) - 1,2) * np.power(m.sin(angle),2) + np.power(m.cos(angle),2),1/2)

print(E)











#semi-major axis
# a = 1/((2/dist) - (velocity*velocity/(G*M)) )
# print(dist,a)
#eccentricity vector



vel1 = -velocity*vel
pos = dist*pos

h = np.cross(pos,vel1)
H = np.cross(vel1,h)

print(pos/dist)
evec = (H/(G*M)) - (pos/dist)
print(np.linalg.norm(evec))


#test1

#true anomaly
# factor = dist*velocity*velocity/(G*M)
# tanv = factor*m.sin(angle)*m.cos(angle)/(factor*np.power(m.sin(angle),2) - 1)
# v = m.atan(tanv)
# print(v*180/np.pi)


# eccentricity
v=(90+abs(DEC[i]))*np.pi/180
u=RA[i]*np.pi/180
x=np.cos(u)*np.sin(v)
y=np.sin(u)*np.sin(v)
z=np.cos(v)


v=(90+abs(DEC[j]))*np.pi/180
u=RA[j]*np.pi/180
X=np.cos(u)*np.sin(v)
Y=np.sin(u)*np.sin(v)
Z=np.cos(v)

velv = (X-x,Y-y,Z-z)                 
posv = (x,y,z)

velvm = np.linalg.norm(velv)
posvm = np.linalg.norm(posv)

vel = velv/velvm
pos = posv/posvm

dt = tl[j] - tl[i]
# cubic root of mass of the earth and gravitational constant
velocity = np.power((np.linalg.norm(velv)/dt)*6.674*np.power(10,-11.0)*5.972*np.power(10,24.0),1/3)
dist = velocity/(np.linalg.norm(velv)/dt)

vel_vec =velocity*vel
pos_vec = dist*pos
mu = G*M

ecc_vec = ((velocity*velocity*pos_vec)/(mu)) -  ((np.dot(pos_vec,vel_vec)*vel_vec)/(mu)) - ((pos_vec)/(dist))
print(np.linalg.norm(ecc_vec))

h_vec = np.cross(pos_vec,vel_vec)
h_mag = np.linalg.norm(h_vec)
i = m.acos(h_vec[2]/h_mag)
print(i*180/np.pi)
k_vec = (0,0,1)
n_vec = np.cross(k_vec,h_vec)
n_mag = np.linalg.norm(n_vec)
ascen_node = m.acos(n_vec[0]/n_mag)*180/np.pi
if(n_vec[1]<0):
    ascen_node = 180 +ascen_node
print(ascen_node)











#ascending node
mag=4
pv1 = (6378*x,6378*y,6378*z)
print(50000*(z-Z)/velvm,6378*z)
if (z >=0):
   
    
    for mag in np.arange(0,10000,1):
        summ = 6378*z + mag*(Z-z)/velvm
        if (abs(summ) <= 1):
            pv3 = (6378*x + mag*(X-x)/velvm, 6378*y + mag*(Y-y)/velvm, 0)
    
    
if (z<=0):

    print('yes')
    for mag in np.arange(0,50000,1):
        summ = 6378*z + mag*(Z-z)/velvm
        
        if (abs(summ) <= 1):
            pv3 = (6378*x + mag*(X-x)/velvm, 6378*y + mag*(Y-y)/velvm, 0)
            
   
pv3m = np.linalg.norm(pv3)
pv3u = pv3/pv3m

#vernal equinox
pv4 = (1,0,0)

angle = m.acos(np.dot(pv3u,pv4))
print('Ascending node  ' + str(180 + (angle*180/np.pi)) + '  degrees')
# 2458950.72211638












#transformation between cartestian to spherical
v=(90+abs(DEC[2]))*np.pi/180
u=RA[2]*np.pi/180
x=np.cos(u)*np.sin(v)
y=np.sin(u)*np.sin(v)
z=np.cos(v)


v=(90+abs(DEC[3]))*np.pi/180
u=RA[3]*np.pi/180
X=np.cos(u)*np.sin(v)
Y=np.sin(u)*np.sin(v)
Z=np.cos(v)



velv = (x-X,y-Y,z-Z)                 
posv = (x,y,z)

velvm = np.linalg.norm(velv)
posvm = np.linalg.norm(posv)

vel = velv/velvm
pos = posv/posvm


vec1 = (X-x,Y-y,Z-z)
dt = tl[3] - tl[2]
# cubic root of mass of the earth and gravitational constant
velocity = np.power((np.linalg.norm(vec1)/dt)*6.674*np.power(10,-20.0)*5.972*np.power(10,24.0),1/3)
dist = velocity/(np.linalg.norm(vec1)/dt)

print(dist)

v=(90)*np.pi/180
u=RA[2]*np.pi/180
X0=dist*np.cos(u)*np.sin(v)
Y0=dist*np.sin(u)*np.sin(v)
Z0=dist*np.cos(v)






v=(90+abs(DEC[30]))*np.pi/180
u=RA[30]*np.pi/180
x=np.cos(u)*np.sin(v)
y=np.sin(u)*np.sin(v)
z=np.cos(v)


v=(90+abs(DEC[31]))*np.pi/180
u=RA[31]*np.pi/180
X=np.cos(u)*np.sin(v)
Y=np.sin(u)*np.sin(v)
Z=np.cos(v)



velv = (x-X,y-Y,z-Z)                 
posv = (x,y,z)

velvm = np.linalg.norm(velv)
posvm = np.linalg.norm(posv)

vel = velv/velvm
pos = posv/posvm


vec1 = (X-x,Y-y,Z-z)
dt = tl[31] - tl[30]
# cubic root of mass of the earth and gravitational constant
velocity = np.power((np.linalg.norm(vec1)/dt)*6.674*np.power(10,-20.0)*5.972*np.power(10,24.0),1/3)
dist = velocity/(np.linalg.norm(vec1)/dt)
print(dist)

v=(90)*np.pi/180
u=RA[30]*np.pi/180
X1=dist*np.cos(u)*np.sin(v)
Y1=dist*np.sin(u)*np.sin(v)
Z1=dist*np.cos(v)


print(X1,Y1)
# X1=-3519.333454955899
# Y1 = 5820.427417263025
# X0 = 624.74125988126
# Y0 = -6789.5595257

#ellipse equation
m1 = abs((Y0*Y0*X1*X1 - X0*X0*Y1*Y1)/(X1*X1 - X0*X0))

b = np.power(m1,0.5)
m2 =abs((Y1*Y1*X0*X0 - X1*X1*Y0*Y0)/(Y1*Y1 - Y0*Y0))

a = np.power(m2,0.5)


print(a,b)


# E=np.power((a*a- b*b),0.5)/a
# print(E)



























