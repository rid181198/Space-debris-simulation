# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 11:48:25 2020

@author: Admin
"""


import matplotlib.pyplot as plt
import numpy as np
import config as c
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
import math as m
import Noises as ns
import qe as QE

line1 = c.line_1
line2 = c.line_2
satellite = twoline2rv(line1, line2, wgs84)
position, velocity=satellite.propagate(c.data[0],c.data[1],c.data[2],c.data[3],c.data[4],c.data[5])

class earth():
        def __init__(self,rx,ry,rz):
            self.rx=rx
            self.ry=ry
            self.rz=rz
        
        #earth's parametric equation
        def eaplot(self):
            #array of input variables
            u=np.linspace(0,2*np.pi,180)
            v=np.linspace(0,np.pi,180)
            
            #equation of ellipsoid with different radius 
            x=self.rx*np.outer(np.cos(u),np.sin(v))
            y=self.ry*np.outer(np.sin(u),np.sin(v))
            z=self.rz*np.outer(np.ones_like(u),np.cos(v))
            # ax.plot_surface(x,y,z,color='w',rstride=4,cstride=4,alpha=0.2)

class vector():
        array=[]
        def __init__(self,x,y,z):
            self.x=x
            self.y=y
            self.z=z
        
        #vector defination
        def plots(self):
            vector.array = [self.x,self.y,self.z]
            t=np.arange(0,1,0.1)
            X = self.x * t
            Y = self.y * t
            Z = self.z * t
            # ax.plot(X,Y,Z,color='b')
earth_wgs84 = earth(6378.137,6356.752314245,6378.137)

class lat_lon():
        vector_ll=[]
        def __init__(self,lat,lon,theta):
            self.lat=lat
            self.lon=lon
            self.theta=theta
        
        #parametric equation for lat_lon
        def location(self):
            u=np.pi*(self.lon)/180
            v=np.pi*(90-self.lat)/180
    
        #equation of ellipsoid with different radius 
            x=earth_wgs84.rx*np.cos(u)*np.sin(v)
            y=earth_wgs84.ry*np.sin(u)*np.sin(v)
            z=earth_wgs84.rz*np.cos(v)
            
            #rotation
            px = x * np.cos(self.theta*np.pi/180) - y * np.sin(self.theta*np.pi/180)
            py = x * np.sin(self.theta*np.pi/180) + y * np.cos(self.theta*np.pi/180)
            pz = z
            lat_lon.vector_ll=[px,py,pz]
            # vector_class = vector(px,py,pz)
            # vector_class.plots()

class unit_vector():
        array=[]
        def __init__(self,x,y,z):
            self.x=x
            self.y=y
            self.z=z
        def convert(self):
            vector.array = [self.x,self.y,self.z]
            magnitude = np.linalg.norm(vector.array)
            unit = vector.array/magnitude
            return unit            

#axis directed towards vernal equinox
vernal_vector = vector(earth_wgs84.rx,0,0)
    
#rotation axis of earth
zaxis_vector = vector(0,0,earth_wgs84.rz)
    
#poistion vector of satellite
sat_vector = vector(position[0],position[1],position[2])


#static declination and right ascension
class ra_dec(): 
        #time formate yyyy,mm,dd,HH,MM,SS  
        declination =0 
        RA =0
        satposition=[]
        def position(*args):
       
            line1 =c.line_1
            line2 = c.line_2
            satellite = twoline2rv(line1, line2, wgs84)
            position, velocity=satellite.propagate(args[0],args[1],args[2],args[3],args[4],args[5])
            ra_dec.satposition = (position[0],position[1],position[2])
            #poistion vector of satellite
            #sat_vector = vector(position[0],position[1],position[2])
            sat_unit = unit_vector(position[0],position[1],position[2]).convert()
    
            #satellite position projection on equator plane
            #sat_xyvector = vector(position[0],position[1],0)
            sat_xyunit = unit_vector(position[0],position[1],0).convert()
            
            if (position[2]>=0):    
                declination = m.acos(np.dot(sat_unit, sat_xyunit))*180/np.pi
            if (position[2]<0):
                declination = -m.acos(np.dot(sat_unit, sat_xyunit))*180/np.pi
            ra_dec.declination = declination
            
            
            #RA
            vernal_unit = unit_vector(earth_wgs84.rx,0,0).convert()
            RA =  360 -((m.acos(np.dot(sat_xyunit,vernal_unit)))*180/np.pi)   
            if (position[0]<0 and position[1]>0):
                RA = (((m.acos(np.dot(sat_xyunit,vernal_unit)))*180/np.pi))
            if (position[0]>0 and position[1]>0):
                RA = ((m.acos(np.dot(sat_xyunit,vernal_unit)))*180/np.pi)
            ra_dec.RA=RA


#static altitude azimuth
class alt_az():
    
         #time formate yyyy,mm,dd,HH,MM,SS, Lat, Lon
         altitude = 0
         azimuth =0
         JD=0
         elevation=0
         LST=0
         location_obs = []
         satposition = []
         def position(*args):
             ra_dec.position(args[0],args[1],args[2],args[3],args[4],args[5])
             RA = ra_dec.RA
             Dec = ra_dec.declination
             yyyy = args[0]
             mm = args[1]
             dd = args[2]
             HH =args[3]
             MM = args[4]
             SS = args[5]
            
             JD = (367*yyyy - int((7*(yyyy + int((mm + 9 )/12))/4 ))- int(3*(int((yyyy+(mm-9)/7 )/100)+1)/4)   +   int(275*mm/9) +dd + 1721028.5 ) + (HH+(MM/60) + (SS/3600) )/24
             alt_az.JD =JD
             Lat = args[6]
             Long = args[7]
             #time for J2000
             day = (JD - 2451545) 
             LST = (100.46 + 0.985647 * day + Long + 15 * (HH + MM / 60) + 360) - (((int)((100.46 + 0.985647 * day + Long + 15 * (HH + MM / 60) + 360)/360))*360)
             alt_az.LST=LST
             HA = (LST - RA + 360)- ((int)((LST - RA + 360)/360))*360 
            
             x = m.cos(HA * (np.pi / 180)) * m.cos(Dec * (np.pi / 180))
             y = m.sin(HA * (np.pi / 180)) *m.cos(Dec * (np.pi / 180))
             z = m.sin(Dec * (np.pi / 180))
              
             xhor = x * m.cos((90 - Lat) * (np.pi / 180)) - z*m.sin((90 - Lat) * (np.pi / 180))
             yhor = y
             zhor = x * m.sin((90 - Lat) * (np.pi / 180)) + z*m.cos((90 - Lat) * (np.pi / 180))
        
             if (xhor>=0 and yhor>=0):
                 azimuth= 180 +((m.atan((yhor/xhor))) * (180 / np.pi)) 
             if (xhor<0 and yhor>0):
                 azimuth=360+((m.atan((yhor/xhor))) * (180 / np.pi)) 
             if (xhor<0 and yhor<0):
                 azimuth=((m.atan((yhor/xhor))) * (180 / np.pi)) 
             if (xhor>0 and yhor<0):
                 azimuth=180+((m.atan((yhor/xhor))) * (180 / np.pi)) 
             altitude = (m.asin(zhor) )* (180 / np.pi)
             alt_az.altitude = altitude
             alt_az.azimuth = azimuth
             
             
             #elevation
             #satellite oberver relative position
             alt_az.location_obs = lat_lon(c.data[6],c.data[7],LST-Long)
             alt_az.location_obs.location()
             sat_obsunit = unit_vector(ra_dec.satposition[0]-alt_az.location_obs.vector_ll[0],ra_dec.satposition[1]-alt_az.location_obs.vector_ll[1],ra_dec.satposition[2]-alt_az.location_obs.vector_ll[2]).convert()
             obs_loc = unit_vector(alt_az.location_obs.vector_ll[0],alt_az.location_obs.vector_ll[1],alt_az.location_obs.vector_ll[2]).convert()
             elevation = 90-(m.acos(np.dot(obs_loc,sat_obsunit))*180/np.pi)
             alt_az.elevation = elevation
    
             alt_az.satposition = (ra_dec.satposition[0],ra_dec.satposition[1],ra_dec.satposition[2])


class photons():
    adc =0
    num=0
    d=0
    #time s, size m, reflectivity,apparent mag.,exposure in s
    def get(*args):
            alt_az.position(c.data[0],c.data[1],c.data[2],c.data[3],c.data[4],c.data[5]+args[0],c.data[6],c.data[7])
                           
                          
                               
            sun_radius = 149597870
            yyyy = c.data[0]
            mm = c.data[1]
            dd = c.data[2]
            HH =c.data[3]
            MM = c.data[4]
            SS = c.data[5]+args[0]
            
            JD = (367*yyyy - int((7*(yyyy + int((mm + 9 )/12))/4 ))- int(3*(int((yyyy+(mm-9)/7 )/100)+1)/4)   +   int(275*mm/9) +dd + 1721028.5 ) + (HH+(MM/60) + (SS/3600) )/24
            
            #time for J2000
            n = (JD - 2451545) 
            L = 280.460 + 0.9856474*n
            g = 357.528 + 0.9856003*n
            e = (23.439 - 0.0000004*n)
            g = g*np.pi/180
            lam = (L + 1.915*np.sin(g) + 0.020*np.sin(2*g))
            
            lam = lam*np.pi/180
            R = sun_radius*(1.00014 - 0.01671*np.cos(g) - 0.00014*np.cos(2*g))
            e = e*np.pi/180
            X= R*(np.cos(lam))
            Y=R*(np.cos(e)*np.sin(lam))
            Z=R*(np.sin(e)*np.sin(lam))
            sun_v = [X,Y,Z]
            sat_v = (alt_az.satposition[0],alt_az.satposition[1],alt_az.satposition[2])
            
            obs_v= (alt_az.location_obs.vector_ll[0],alt_az.location_obs.vector_ll[1],alt_az.location_obs.vector_ll[2])
            
            sat_obs = (alt_az.satposition[0]-alt_az.location_obs.vector_ll[0],alt_az.satposition[1]-alt_az.location_obs.vector_ll[1],alt_az.satposition[2]-alt_az.location_obs.vector_ll[2])
            sat_sun = (alt_az.satposition[0]-X,alt_az.satposition[1]-Y,alt_az.satposition[2]-Z)
            sun_obs = (X-alt_az.location_obs.vector_ll[0],Y-alt_az.location_obs.vector_ll[1],Z-alt_az.location_obs.vector_ll[2])
            sat_obs_m = np.linalg.norm(sat_obs)
            sat_sun_m = np.linalg.norm(sat_sun)
            sun_obs_m = np.linalg.norm(sun_obs)
            
            phi = (m.acos(((np.power(sat_obs_m,2))+(np.power(sat_sun_m,2)) - (np.power(sun_obs_m,2)))/(2*sat_obs_m*sat_sun_m) ))
            # phi = np.pi - phi
            
            F = (2/(3*np.pi*np.pi) )*((np.pi-phi)*np.cos(phi) + np.sin(phi))
            
            d=(args[1]/1000)
            A = (np.pi*d*d)/4
            r=args[2]
            
            D=sat_obs_m
                               
            elevation = alt_az.elevation
                               
                               
            #apparent magnitude
            ze = (90-elevation)*np.pi/180
            
            am = np.power(np.cos(ze),-1) #- 0.0018167*(np.power(np.cos(ze),-1) - 1) - 0.002875*np.power((np.power(np.cos(ze),-1) - 1),2) - 0.0008083*np.power((np.power(np.cos(ze),-1) - 1),3)
            
            mo =-26.8
            
            
            #atmospheric extinction
            h=0.018
            wav = 5440/10000
            kr = ((0.0095*m.exp(-h/8))/np.power(wav,4))*(0.23465 + (107.6/(146-np.power(wav,-2))) + (0.93161/(41-np.power(wav,-2))))
            ka = 0.087*m.exp(-h/1.5)/(np.power(wav,0.8))
            ko = (839.4375*m.exp(-131*(wav-0.26))) + 0.03811562*m.exp(-188*(wav-0.59)*(wav-0.59))
            
            K = kr + ka + ko
            
            mv= float(mo -2.5*(m.log(A*r*F/(D*D),10)))
            mv = mv+ K*am
                     
            # r = np.power(10,(18+26.58)/(-2.5))*D*D/(A*F)
            # print(mv)
            
            #absolute magnitude
            M =mv-(5*(m.log(D)))
            D = 3.24078e-14*D
            M =mv-(5*(m.log(D,10)))+5
            # print(M)
           
    
            #photon numbers using sun
            # fx = 1373*np.power(10,(-26.8-mv)/2.5)
            
            #photon numbers using vega star in erg/cm^2*s*A
            fx = (363.1)*np.power(10,-11.0)*np.power(10,-mv/2.5)
            # fx = (2.180100)*np.power(10,-9.0)*np.power(10,-mv/2.5)*np.power(10,-4.0)
            
            #wavelength in Angstrom.
            wavelength_band = 880
            wavelength_central = 5440
            
            #integration time t in sec.
            it=args[3]
            
            #area of telescope in cm^2
            area = np.pi*(1.2/2)*(1.2/2)*10000
            #pixel area
            # area = 9*9*(np.power(10,-8.0))
            
            num = fx*wavelength_central*wavelength_band*area*it*(np.power(10,-7.0))*(np.power(10,-10.0))/(299792458*6.626176*np.power(10,-34.0))
            
            #quantum efficiency
            for k in QE.qe:
                 diff_wav = (wavelength_central/10) - k
                 if (0<= diff_wav <=10):
                     qeff = QE.qe[k]
                                       
            photons.num = qeff*num/100
           
            
                                                           
            
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
            fov2 = no_pixel*pp*0.000277778
            pixel = 1024
            P = fov2/pixel
            #size of the satellite in km
            size = d
            #distance in km
            photons.d =sat_obs_m
            #angle
            alpha = size*180/(photons.d*np.pi)
            pcoverage = alpha/P
            
            
            
            #noise
            #argument as follows pixels, exposure time, quantum efficiency, area of pixel in micron, photon flux      
            ns.noise.get(pcoverage*4.5,it,50,169,photons.num)
            noise = ns.noise.NT
       
            
            #brightness
            #pogsun formula
            # m=m2 -2.5*np.log(B1/B2,10)
            #m2, B2 are reference star's magnitude and brightness
            
            #rough ADC converter
            photons.adc = int(65355*(photons.num-noise)/(100000*pcoverage*4.5*700))
            # photons.adc = int(16777215*(photons.num-noise)/(1000000*pcoverage*4.5))        
        
    
    
ra_dec.position(c.data[0],c.data[1],c.data[2],c.data[3],c.data[4],c.data[5])
alt_az.position(c.data[0],c.data[1],c.data[2],c.data[3],c.data[4],c.data[5],c.data[6],c.data[7])

sat_vector = (position[0],position[1],position[2])
obs_vector = (alt_az.location_obs.vector_ll[0],alt_az.location_obs.vector_ll[1],alt_az.location_obs.vector_ll[2])
    
sat_obs_vector = (position[0]-alt_az.location_obs.vector_ll[0],position[1]-alt_az.location_obs.vector_ll[1],position[2]-alt_az.location_obs.vector_ll[2])
distance = np.linalg.norm(sat_obs_vector)







