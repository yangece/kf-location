# -*- coding: utf-8 -*-
# input from subspace method

from math import *
# from dfl import *
# from spin import *
from pylab import *
from time import sleep
from quickSave import *
import cPickle

# for plot
ion()
params = { 'axes.labelsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18}           
rcParams.update(params)

ftm = 0.3048
pi = 3.1415926

vrti_map_file = './House_sub_yang_k36.txt' # use map
vrti_interp_file = './House_interp_yang.txt'
configFile = './house_may_2010.cfg'
#fig_name = 'track_sub_yang_Q2_Cx1.eps'


configDict = quickLoad(configFile)
nodes = configDict['nodes']
nodes = nodes*ftm
nodeX = nodes[:,0]
nodeY = nodes[:,1]

########################
# Kalman filtering
########################

# read measurement vector from data file
file_s = open(vrti_map_file, 'r')
Sest = cPickle.load(file_s)
file_s.close()
Sest = array(Sest)*ftm;
Xest = Sest[:,0]; Yest = Sest[:,1]

file_i = open(vrti_interp_file, 'r')
interp = cPickle.load(file_i)
file_i.close()
interp = array(interp)*ftm;
interpX = interp[:,0]; interpY = interp[:,1]

# initialization
startP = 15 # remove the first few bad rti estimations: start from number 15
#startP = 20
mu_s = array([Xest[startP],Yest[startP],cos(pi),sin(pi)])

Cs = 100. * eye(4);

dt = 1.; 
#dt = .6;

A = array([[1, 0, dt, 0], [0, 1, 0, dt], [0, 0, 1, 0], [0, 0, 0, 1]]);

#var_u = 100.; var_v = 100.; # larger: better for agile motions
var_u = 2.; var_v = 2.; # smaller: less process noise, better for straight line path


# Q: the variance of the process noise
Q = array([ [0,0,0,0], [0,0,0,0], [0,0,var_u,0], [0,0,0,var_v] ]); 

H = array([ [1,0,0,0], [0,1,0,0] ]);

# variance of the measurement noise
var_x = 5.; var_y = 5. ### measurement noise

# observation variance
Cx = array([ [var_x, 0], [0, var_y] ])

flag_anim = False;
flag_anim = True;
# first guess for s_est
s_est = mu_s;
MSE = Cs;

### start iteration

# S: estimated state vector
S_x=[]; S_y=[]; S_std=[]; S_se=[]; 
I_x=[]; I_y=[]
KF_s = []; KF_i = []; # for KF output file
for i in range(startP,len(Xest)): # from start point = 5
  
    # get measurement vector
    x = array([Xest[i], Yest[i]]);
    
    # prediction
    s_pred = dot(A,s_est);
    
    # minimum prediction MSE
    MSE_p = dot(A, dot(MSE, A.T)) + Q;
    
    # Kalman gain
    K1 = inv(Cx + dot(H, dot(MSE_p, H.T)));
    K = dot(MSE_p, dot(H.T, K1));
    
    # correction
    s_est = s_pred + dot(K, (x - dot(H,s_pred)));
    
    # minimum MSE
    MSE = dot((eye(4) - dot(K,H)), MSE_p);
    
    s_est_p = array([s_est[0],s_est[1]])
    KF_s.append(s_est_p); KF_i.append(interp[i])
    S_x.append(s_est[0]); S_y.append(s_est[1]); 
    I_x.append(interpX[i]); I_y.append(interpY[i])

# plot the KF tracking results

for i in range(0,len(Xest)-startP): # remove startP
    #if i == 0:
        #plot([S_x[i]],[S_y[i]],'ro')
    #else:
        #plot([S_x[i-1],S_x[i]], [S_y[i-1],S_y[i]],'ro')
        #plot([S_x[i-1],S_x[i]], [S_y[i-1],S_y[i]],'r')        
        #plot([I_x[i]],[I_y[i]],'bo')        
    
    s_std = sqrt((I_x[i]-S_x[i])**2 + (I_y[i]-S_y[i])**2)
    S_std.append(s_std)
    #print s_std
    
    s_se = (I_x[i]-S_x[i])**2 + (I_y[i]-S_y[i])**2
    S_se.append(s_se)
    

rmseD = sqrt(mean(S_se))
print 'rmse', rmseD
#print 'mean', mean(S_std)


s_rti = array(KF_s)
rti_x = s_rti[:,0]; rti_y = s_rti[:,1]
s_interp = array(KF_i)
interpX = s_interp[:,0]; interpY = s_interp[:,1]

params = { 'axes.labelsize': 18,
           'xtick.labelsize': 18,
           'ytick.labelsize': 18,
           'legend.fontsize': 15}
           
rcParams.update(params)

figure(2)
plot([.6,8.2],[1.8,1.8],'k-',linewidth=4, label='Walls')
plot([.6,.6],[1.8,7.5],'k-',linewidth=4)
plot([.6,10],[7.5,7.5],'k-',linewidth=4)
plot([8.2,8.2],[1.8,7.5],'k-',linewidth=4)
plot([2.4,2.4],[7.5,9.],'k-',linewidth=4)

plot([8.2,8.2],[5.7,6.6],'k-',linewidth=8)  # right door
plot([3.1,3.8],[7.5,7.5],'k-',linewidth=8)
nodePlot, = plot(nodes[:,0],nodes[:,1],'k.',markersize=20)

knownPlot = plot(interpX[startP:len(interpX)], interpY[startP:len(interpY)],'g:', linewidth=2)
rtiPlot = plot(rti_x[startP: len(rti_x)], rti_y[startP: len(rti_y)],'b-', linewidth=2)
# legend((knownPlot, rtiPlot), ('Known path','Tracking estimates'), numpoints=1, loc=(0.35,0.92))
xlim((-1,10)),ylim((-1,9))
xlabel('X (m)'); ylabel('Y (m)')

savefig("kf.png")


#file_s = open("KF_s_yang_sub_k0_test.txt", 'w')
#cPickle.dump(KF_s, file_s)
#file_s.close()

#file_i = open("KF_i_yang_Cx200_Q400.txt", 'w')
#cPickle.dump(KF_i, file_i)
#file_i.close()
