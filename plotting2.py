

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib import rcParams 
from numpy import linalg as LA

x = np.linspace(0.2, 5, 100) 

errProb1 = np.loadtxt("errProb1.dat")
errRK2 = np.loadtxt("errRK2.dat")
errRK4 = np.loadtxt("errRK4.dat")
errRKf45 = np.loadtxt("errRKf45.dat")












# In[2]:

fig1 = plt.figure()
plt.plot(errProb1[:,0],errProb1[:,1] ,'r-o')
plt.xscale('log')
# to set y-axis to logscale
plt.yscale('log')
plt.xlabel('number of iterations')
plt.ylabel('error, N')
plt.title("error vs iterations, for myEuler vs Analytical, problem 1")
plt.savefig("prob1_err_myEuler.pdf")
plt.close(fig1)



fig2 = plt.figure()
plt.plot(errProb1[:6,0],errProb1[:6,2] ,'r-o')
plt.xscale('log')
# to set y-axis to logscale
plt.yscale('log')
plt.xlabel('number of iterations')
plt.ylabel('error, N')
plt.title("error vs iterations, for GSL vs Analytical, problem 1")
plt.savefig("prob1_err_GSL.pdf")
plt.close(fig2)


fig3 = plt.figure()
plt.plot(errRK2[:5,0],errRK2[:5,1] ,'r-o')
plt.xscale('log')
# to set y-axis to logscale
plt.yscale('log')
plt.xlabel('number of iterations')
plt.ylabel('error, N')
plt.title("error vs iterations, for rk2 vs accurate solution, problem 2")
plt.savefig("prob2_err_rk2.pdf")
plt.close(fig2)

fig4 = plt.figure()
plt.plot(errRK4[:4,0],errRK4[:4,1] ,'r-o')
plt.xscale('log')
# to set y-axis to logscale
plt.yscale('log')
plt.xlabel('number of iterations')
plt.ylabel('error, N')
plt.title("error vs iterations, for rk4 vs accurate solution, problem 2")
plt.savefig("prob2_err_rk4.pdf")
plt.close(fig2)

fig5 = plt.figure()
plt.plot(errRKf45[1:5,0],errRKf45[1:5,1] ,'r-o')
plt.xscale('log')
# to set y-axis to logscale
plt.yscale('log')
plt.xlabel('number of iterations')
plt.ylabel('error, N')
plt.title("error vs iterations, for rk45 vs accurate solution, problem 2")
plt.savefig("prob2_err_rkf45.pdf")
plt.close(fig2)










