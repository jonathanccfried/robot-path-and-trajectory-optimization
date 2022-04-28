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
import cvxpy


sys.path.append('../../toolbox')
from robots_def import *
from lambdacalculator import lambda_calculator
from cvx_ldot_optimizer import cvx_ldot_optimizer
from constant_speed_solver import constant_speed_solver
from one_step_optimizer import one_step_q_optimizer
from one_step_optimizer_v3 import one_step_q_optimizer_v3

col_names = ['X', 'Y', 'Z', 'direction_x', 'direction_y', 'direction_z']
data = read_csv("../../data/from_ge/Curve_in_base_frame2.csv", names=col_names)
curve_x = data['X'].tolist()
curve_y = data['Y'].tolist()
curve_z = data['Z'].tolist()
curve_direction_x = data['direction_x'].tolist()
curve_direction_y = data['direction_y'].tolist()
curve_direction_z = data['direction_z'].tolist()
curve = np.vstack((curve_x, curve_y, curve_z)).T
curve_normal = np.vstack((curve_direction_x, curve_direction_y, curve_direction_z)).T

#print(curve.shape)
#print(curve_normal.shape)
lam = lambda_calculator(np.transpose(curve),50000)
#print(lam.shape)
robot1 = abb6640()

base2_R = np.eye(3)
base2_p = np.zeros(3)
joint_vel_limit = np.radians(np.reshape([110, 90, 90, 150, 120, 235],[6,1]))
joint_acc_limit = np.radians(np.reshape([312, 292, 418, 2407, 1547, 3400],[6,1]))
steps = 500

###decrease curve density to simplify computation
num_per_step = int(len(curve) / steps)
#print(num_per_step)
curve=curve[0:-1:num_per_step]
curve_normal=curve_normal[0:-1:num_per_step]
lam = lam[0:-1:num_per_step]


###find path length
theta = 4.00553063332701 * np.ones(len(curve))

def fwd_all(robot,q_all, base_R=np.eye(3), base_p=np.array([0, 0, 0])):
    pose_p_all = []
    pose_R_all = []
    for q in q_all:
        pose_temp = robot.fwd(q, base_R, base_p)
        pose_p_all.append(pose_temp.p)
        pose_R_all.append(pose_temp.R)

    return Transform_all(pose_p_all, pose_R_all)

def direction2R(v_norm, v_tang):
    v_norm = v_norm / np.linalg.norm(v_norm)
    theta1 = np.arccos(np.dot(np.array([0, 0, 1]), v_norm))
    ###rotation to align z axis with curve normal
    axis_temp = np.cross(np.array([0, 0, 1]), v_norm)
    R1 = rot(axis_temp / np.linalg.norm(axis_temp), theta1)

    ###find correct x direction
    v_temp = v_tang - v_norm * np.dot(v_tang, v_norm) / np.linalg.norm(v_norm)

    ###get as ngle to rotate
    theta2 = np.arccos(np.dot(R1[:, 0], v_temp / np.linalg.norm(v_temp)))

    axis_temp = np.cross(R1[:, 0], v_temp)
    axis_temp = axis_temp / np.linalg.norm(axis_temp)

    ###rotation about z axis to minimize x direction error
    R2 = rot(np.array([0, 0, np.sign(np.dot(axis_temp, v_norm))]), theta2)

    return np.dot(R1, R2)




for i in range(len(curve)):
    if i == 0:
        R_temp = direction2R(curve_normal[i], -curve[i + 1] + curve[i])
        R = np.dot(R_temp, Rz(theta[i]))
        try:
            q_out = [robot1.inv(curve[i], R)[0]]
        except:
            traceback.print_exc()


    else:
        R_temp = direction2R(curve_normal[i], -curve[i] + curve[i - 1])
        R = np.dot(R_temp, Rz(theta[i]))
        try:
            ###get closet config to previous one
            q_inv_all = robot1.inv(curve[i], R)
            temp_q = q_inv_all - q_out[-1]
            order = np.argsort(np.linalg.norm(temp_q, axis=1))
            q_out.append(q_inv_all[order[0]])
        except:
            traceback.print_exc()

q_out = np.asarray(np.transpose(q_out))
f = fwd_all(robot1,np.transpose(q_out))
#print(q_out.shape)
# best = scipy.io.loadmat('best.mat')
# q_out = best['qbest']
# print(q_out.shape)
# lam = best['l']
# lam = lam.reshape(lam.shape[1])
# print(lam.shape)
# qprimemat = best['qprimeorg']
# qdoubleprimemat = best['qdoubleprimeorg']
size = q_out.shape
#print(size[0])
opt_steps = 1
step_size = 0.05
discount = 0.1
n = 15 #polynomial degree
a = np.zeros([size[0],n+1])
da = np.zeros([size[0],n])
dda = np.zeros([size[0],n-1])
qprime = np.zeros([size[0],len(lam)])
qdoubleprime = np.zeros([size[0],len(lam)])
q_fit = q_out
for i in range(size[0]):
    a[i,:] = np.polyfit(lam[:],q_out[i,:],n)
    q_fit[i,] = np.polyval(a[i,:],lam)
    da[i,:] = np.polyder(a[i,:])
    qprime[i,:] = np.polyval(da[i,:],lam)
    dda[i,:] = np.polyder(da[i,:])
    qdoubleprime[i,:] = np.polyval(dda[i,:],lam)
ffit = fwd_all(robot1, np.transpose(q_fit))
#print(qprime.shape)
#print(qdoubleprime.shape)
#print(lam.shape)
# print(a[0,:])
# print(qdoubleprime[0,-1])
# print(q_out[0,-1])
# plt.plot(a[0,:])
# plt.show()

poli = np.zeros([size[1],n+1])
poliprime = np.zeros([size[1],n])
polidoubleprime = np.zeros([size[1],n-1])
for i in range(size[1]):
    for j in range(n+1):
        poli[i,j] = lam[i]**(n-j)
        if j < n:
            poliprime[i,j] = (n-j)*(lam[i]**(n-j-1))
        if j < (n-1):
            polidoubleprime[i,j] = (n-j)*(n-j-1)*(lam[i]**(n-j-2))

polimatrix = np.zeros([size[1]*size[0],size[0]*(n+1)])
for i in range(size[1]):
    for j in range(size[0]):
        polimatrix[size[0]*i+j,j*(n+1):(j+1)*(n+1)] = poli[i,:]

ldotmin, t_final, indexmin, locktype, lag = constant_speed_solver(qprime,qdoubleprime,-joint_vel_limit,joint_vel_limit,-joint_acc_limit,joint_acc_limit)

print('ldotmin=',ldotmin)
print('t_final=',t_final)
print('indexmin=',indexmin)
print('locktype=',locktype)
print('lag=',lag)
q_out2 = q_out



start_time = time.time()
for s in range(opt_steps):
    #cvx_ldot_optimizer(qprime,qdoubleprime,lam,-joint_vel_limit,joint_vel_limit,-joint_acc_limit,joint_acc_limit)
    ldotmin, t_final, indexmin, locktype, lag = constant_speed_solver(qprime,qdoubleprime,-joint_vel_limit,joint_vel_limit,-joint_acc_limit,joint_acc_limit)
    # dtheta = one_step_q_optimizer(robot1,curve_normal,q_out2,lam,a,da,dda,polimatrix,poliprime,polidoubleprime,ldotmin,lag,indexmin,locktype,n,size,taskdof=5,eps=0.5)
    dq = one_step_q_optimizer_v3(robot1, curve_normal, q_out2, lam, a, da, dda, polimatrix, poliprime, polidoubleprime,
                                  ldotmin, lag, indexmin, locktype, n, size, taskdof=5, eps=0.5)
    # dq = one_step_q_optimizer(robot1, curve_normal, q_out2, lam, a, da, dda, polimatrix, poliprime, polidoubleprime,
    #                               ldotmin, lag, indexmin, locktype, n, size, taskdof=5, eps=5)
    q_temp = q_out2
    for i in range(size[0]):
        #a[i, :] = a[i,:] - 0.005*np.transpose(dtheta[i*(n+1):(i+1)*(n+1)])
        # a[i, :] = a[i, :] - 0.005 * np.transpose(dtheta[i, :])
        #q_out2[i,:] = np.polyval(a[i,:],lam)
        q_temp[i, :] = q_out2[i, :] - step_size*dq[i,:]
        a[i, :] = np.polyfit(lam[:], q_temp[i, :], n)
        da[i, :] = np.polyder(a[i, :])
        qprime[i, :] = np.polyval(da[i, :], lam)
        dda[i, :] = np.polyder(da[i, :])
        qdoubleprime[i, :] = np.polyval(dda[i, :], lam)
    ldottemp, t_final_temp, indexmin_temp, locktype_temp, lag_temp = constant_speed_solver(qprime, qdoubleprime, -joint_vel_limit,
                                                                      joint_vel_limit, -joint_acc_limit,
                                                                      joint_acc_limit)
    count = 1
    print(count)
    while ldottemp <= ldotmin:
        for i in range(size[0]):
            # a[i, :] = a[i,:] - 0.005*np.transpose(dtheta[i*(n+1):(i+1)*(n+1)])
            # a[i, :] = a[i, :] - 0.005 * np.transpose(dtheta[i, :])
            # q_out2[i,:] = np.polyval(a[i,:],lam)
            q_temp[i, :] = q_out2[i, :] - step_size*(discount**count) * dq[i, :]
            print(discount**count)
            print(step_size*(discount**count))

            a[i, :] = np.polyfit(lam[:], q_temp[i, :], n)
            da[i, :] = np.polyder(a[i, :])
            qprime[i, :] = np.polyval(da[i, :], lam)
            dda[i, :] = np.polyder(da[i, :])
            qdoubleprime[i, :] = np.polyval(dda[i, :], lam)
        ldottemp, t_final_temp, indexmin_temp, locktype_temp, lag_temp = constant_speed_solver(qprime, qdoubleprime,
                                                                                               -joint_vel_limit,
                                                                                               joint_vel_limit,
                                                                                               -joint_acc_limit,
                                                                                               joint_acc_limit)
        print(ldottemp)
        count = count +1
        #print(count)
        if count ==5:
            import sys

            sys.exit()
    q_out2 = q_temp


print("--- %s seconds ---" % (time.time() - start_time))
f2 = fwd_all(robot1,np.transpose(q_out2))

ldotmin, t_final, indexmin, locktype, lag = constant_speed_solver(qprime,qdoubleprime,-joint_vel_limit,joint_vel_limit,-joint_acc_limit,joint_acc_limit)
print('ldotmin=',ldotmin)
print('t_final=',t_final)
print('indexmin=',indexmin)
print('locktype=',locktype)
print('lag=',lag)
plt.plot(f2.p_all[:,0])
plt.plot(ffit.p_all[:,0])
plt.plot(f2.p_all[:,1])
plt.plot(ffit.p_all[:,1])
plt.plot(f2.p_all[:,2])
plt.plot(ffit.p_all[:,2])
plt.show()
#
xnorm = np.linalg.norm(f2.p_all-ffit.p_all)
print(xnorm)
print(np.amax(np.sqrt(np.square(f2.p_all-ffit.p_all)),axis=0))
plt.plot(lam,(np.linalg.norm(f2.p_all-ffit.p_all,axis=1)))
plt.gca().legend(('e_norm'))
plt.show()

angle_error = np.zeros([size[1]])
print(np.linalg.norm(curve_normal[1][:]))
print(robot1.fwd(q_out[:,0]).R[:,[-1]])
import sys
sys.exit()
for i in range(size[1]):
    # print(q_out[i,:].shape)
    angle_error[i] = np.arccos(curve_normal[i][:]@robot1.fwd(q_out[:,i]).R[:,[-1]])

plt.plot(lam,angle_error)
plt.gca().legend(('e_theta'))
plt.show()



#scipy.io.savemat('qcurveopt5000.mat', {'qcurveopt5000':q_out2})
