# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 19:34:59 2022

@author: j2511
"""

import numpy as np
from pandas import *
import sys, traceback, time
from general_robotics_toolbox import *
import matplotlib.pyplot as plt
import scipy.io
import time


sys.path.append('../../toolbox')
from robots_def import *

col_names=['X', 'Y', 'Z','direction_x','direction_y','direction_z'] 
data = read_csv("../../data/from_ge/Curve_in_base_frame2.csv", names=col_names)
curve_x=data['X'].tolist()
curve_y=data['Y'].tolist()
curve_z=data['Z'].tolist()
curve_direction_x=data['direction_x'].tolist()
curve_direction_y=data['direction_y'].tolist()
curve_direction_z=data['direction_z'].tolist()
curve=np.vstack((curve_x, curve_y, curve_z)).T
curve_normal=np.vstack((curve_direction_x, curve_direction_y, curve_direction_z)).T

print(curve.shape)
print(curve_normal.shape)

robot1 = robot1=abb6640()

base2_R=np.eye(3)
base2_p=np.zeros(3)
joint_vel_limit=np.radians([110,90,90,150,120,235])
steps=50000
  
  		###decrease curve density to simplify computation
num_per_step=int(len(curve)/steps)
print(num_per_step)
#curve=curve[0:-1:num_per_step]
#curve_normal=curve_normal[0:-1:num_per_step]
  
  
  		###find path length
theta = 2.45044*np.ones(len(curve))

def direction2R(v_norm,v_tang):
  		v_norm=v_norm/np.linalg.norm(v_norm)
  		theta1 = np.arccos(np.dot(np.array([0,0,1]),v_norm))
  		###rotation to align z axis with curve normal
  		axis_temp=np.cross(np.array([0,0,1]),v_norm)
  		R1=rot(axis_temp/np.linalg.norm(axis_temp),theta1)
  
  		###find correct x direction
  		v_temp=v_tang-v_norm * np.dot(v_tang, v_norm) / np.linalg.norm(v_norm)
  
  		###get as ngle to rotate
  		theta2 = np.arccos(np.dot(R1[:,0],v_temp/np.linalg.norm(v_temp)))
  
  
  		axis_temp=np.cross(R1[:,0],v_temp)
  		axis_temp=axis_temp/np.linalg.norm(axis_temp)
  
  		###rotation about z axis to minimize x direction error
  		R2=rot(np.array([0,0,np.sign(np.dot(axis_temp,v_norm))]),theta2)
  
  		return np.dot(R1,R2)
start_time = time.time()

for i in range(len(curve)):
    if i==0:
        R_temp=direction2R(curve_normal[i],-curve[i+1]+curve[i])
        R=np.dot(R_temp,Rz(theta[i]))
        try:
            q_out=[robot1.inv(curve[i],R)[0]]
        except:
            traceback.print_exc()
            #return 999

    else:
        R_tempemp=direction2R(curve_normal[i],-curve[i]+curve[i-1])
        R=np.dot(R_temp,Rz(theta[i]))
        try:
 			###get closet config to previous one
            q_inv_all=robot1.inv(curve[i],R)
            temp_q=q_inv_all-q_out[-1]
            order=np.argsort(np.linalg.norm(temp_q,axis=1))
            q_out.append(q_inv_all[order[0]])
        except:
            traceback.print_exc()
            #return 999
    #q_out = np.reshape(q_out, [len(curve),6])
    #print(q_out.shape)
    
    
print("--- %s seconds ---" % (time.time() - start_time))
#print(qcurveshist.shape)
#print(thetahist.shape)
scipy.io.savemat('qcurvelong.mat', {'qcurve':q_out})
#print(q_out)
#print(opt)