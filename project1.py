'''
Alexander Bart and Julie Bosley
Feb 19 2018
PHGN498A - Project 1

N x N Coupled Pendula
'''

import numpy as np
import matplotlib.pyplot as plt
import math
import random as rand
import matplotlib.animation as animation
import time as t
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3

plt.close()

N = 2
steps = 10000
finaltime = 5

k = 10
m = 1
R = 5
g = 98.1

def coupled_pendula_3D(steps,timef,N):

    finaltime = timef
    dt = finaltime/steps
    time = np.linspace(0, finaltime, steps+1)

    theta = np.zeros((N+2,N+2,steps+1))
    phi = np.zeros((N+2,N+2,steps+1))

    thdot = np.zeros((N+2,N+2,steps+1))
    phdot = np.zeros((N+2,N+2,steps+1))

    for i in range(1,N+1):
        for j in range(1,N+1):
            theta[i,j,0] = np.pi/4

    phdot[1,1,0] = 3

    #thdot[1,1,0] = -1
    # phdot[0,2,0] = -1


    I = m*(R**2)

    for tt in range(steps):
        for i in range(1,N+1):
            for j in range(1,N+1):
                if 0 <= theta[i,j,tt]-0 < 0.01: theta[i,j,tt] = -0.01
                if 0 > theta[i,j,tt]-0 > -0.01: theta[i,j,tt] = 0.01
                #print(i,j)
                ks = [k,k,k,k] # One for each possible spring per pendulum
                if j+1 > N: ks[0] = 0
                if j-1 < 1: ks[1] = 0
                if i+1 > N: ks[2] = 0
                if i-1 < 1: ks[3] = 0
                #print(i,j,ks)

                thdot[i,j,tt+1] = thdot[i,j,tt] + dt*(phdot[i,j,tt]**2 \
                *np.sin(theta[i,j,tt])*np.cos(theta[i,j,tt]) - \
                g/R*np.sin(theta[i,j,tt]) - 1/m*( \

                ks[0]*(np.sin(theta[i,j+1,tt])*np.cos(phi[i,j+1,tt]) - \
                np.sin(theta[i,j,tt])*np.cos(phi[i,j,tt]))* \
                (-np.cos(theta[i,j,tt])*np.cos(phi[i,j,tt])) + \

                ks[0]*(np.sin(theta[i,j+1,tt])*np.sin(phi[i,j+1,tt]) - \
                np.sin(theta[i,j,tt])*np.sin(phi[i,j,tt]))* \
                (-np.cos(theta[i,j,tt])*np.sin(phi[i,j,tt])) + \

                ks[1]*(np.sin(theta[i,j,tt])*np.cos(phi[i,j,tt]) - \
                np.sin(theta[i,j-1,tt])*np.cos(phi[i,j-1,tt]))*\
                (np.cos(theta[i,j,tt])*np.cos(phi[i,j,tt])) + \

                ks[1]*(np.sin(theta[i,j,tt])*np.sin(phi[i,j,tt]) - \
                np.sin(theta[i,j-1,tt])*np.sin(phi[i,j-1,tt]))*\
                (np.cos(theta[i,j,tt])*np.sin(phi[i,j,tt])) + \

                ks[2]*(np.sin(theta[i+1,j,tt])*np.sin(phi[i+1,j,tt]) - \
                np.sin(theta[i,j,tt])*np.sin(phi[i,j,tt]))*\
                (-np.cos(theta[i,j,tt])*np.sin(phi[i,j,tt])) + \

                ks[2]*(np.sin(theta[i+1,j,tt])*np.cos(phi[i+1,j,tt]) - \
                np.sin(theta[i,j,tt])*np.cos(phi[i,j,tt]))*\
                (-np.cos(theta[i,j,tt])*np.cos(phi[i,j,tt])) + \

                ks[3]*(np.sin(theta[i,j,tt])*np.sin(phi[i,j,tt]) - \
                np.sin(theta[i-1,j,tt])*np.sin(phi[i-1,j,tt]))*\
                (np.cos(theta[i,j,tt])*np.sin(phi[i,j,tt])) + \

                ks[3]*(np.sin(theta[i,j,tt])*np.cos(phi[i,j,tt]) - \
                np.sin(theta[i-1,j,tt])*np.cos(phi[i-1,j,tt]))*\
                (np.cos(theta[i,j,tt])*np.cos(phi[i,j,tt])) ))

                theta[i,j,tt+1] = dt*thdot[i,j,tt] + theta[i,j,tt]


                phdot[i,j,tt+1] = phdot[i,j,tt] + \
                dt*(1/I/(np.sin(theta[i,j,tt])) * (-2*I*phdot[i,j,tt]*\
                thdot[i,j,tt] - k*(R)*( \

                ks[0]*(np.sin(theta[i,j+1,tt])*np.cos(phi[i,j+1,tt]) - \
                np.sin(theta[i,j,tt])*np.cos(phi[i,j,tt]))* \
                (np.sin(phi[i,j,tt])) + \

                ks[0]*(np.sin(theta[i,j+1,tt])*np.sin(phi[i,j+1,tt]) - \
                np.sin(theta[i,j,tt])*np.sin(phi[i,j,tt]))* \
                (-np.cos(phi[i,j,tt])) + \

                ks[1]*(np.sin(theta[i,j,tt])*np.cos(phi[i,j,tt]) - \
                np.sin(theta[i,j-1,tt])*np.cos(phi[i,j-1,tt]))*\
                (-np.sin(phi[i,j,tt])) + \

                ks[1]*(np.sin(theta[i,j,tt])*np.sin(phi[i,j,tt]) - \
                np.sin(theta[i,j-1,tt])*np.sin(phi[i,j-1,tt]))*\
                (np.cos(phi[i,j,tt])) + \

                ks[2]*(np.sin(theta[i+1,j,tt])*np.sin(phi[i+1,j,tt]) - \
                np.sin(theta[i,j,tt])*np.sin(phi[i,j,tt]))*\
                (-np.cos(phi[i,j,tt])) + \

                ks[2]*(np.sin(theta[i+1,j,tt])*np.cos(phi[i+1,j,tt]) - \
                np.sin(theta[i,j,tt])*np.cos(phi[i,j,tt]))*\
                (np.sin(phi[i,j,tt])) + \

                ks[3]*(np.sin(theta[i,j,tt])*np.sin(phi[i,j,tt]) - \
                np.sin(theta[i-1,j,tt])*np.sin(phi[i-1,j,tt]))*\
                (np.cos(phi[i,j,tt]))  + \

                ks[3]*(np.sin(theta[i,j,tt])*np.cos(phi[i,j,tt]) - \
                np.sin(theta[i-1,j,tt])*np.cos(phi[i-1,j,tt]))*\
                (-np.sin(phi[i,j,tt]))  \
                )))

                phi[i,j,tt+1] = dt*phdot[i,j,tt] + phi[i,j,tt]

    energy = []
    for tt in range(steps):
        T = 0
        U = 0
        for i in range(1,N+1):
            for j in range(1,N+1):
                T += (I/2 * thdot[i,j,tt]**2 + \
                0.5*m*(R*np.sin(theta[i,j,tt]))**2 * phdot[i,j,tt]**2)

                ks = [k,k] # One for each possible spring per pendulum
                if j+1 > N: ks[0] = 0
                if i+1 > N: ks[1] = 0

                U += (m*g*R*(1-np.cos(theta[i,j,tt])) + \
                ks[0]/2*(R*np.sin(theta[i,j+1,tt])*np.cos(phi[i,j+1,tt]) - \
                R*np.sin(theta[i,j,tt])*np.cos(phi[i,j,tt]))**2 + \

                ks[0]/2*(R*np.sin(theta[i,j+1,tt])*np.sin(phi[i,j+1,tt]) - \
                R*np.sin(theta[i,j,tt])*np.sin(phi[i,j,tt]))**2 + \

                ks[1]/2*(R*np.sin(theta[i+1,j,tt])*np.cos(phi[i+1,j,tt]) - \
                R*np.sin(theta[i,j,tt])*np.cos(phi[i,j,tt]))**2 + \

                ks[1]/2*(R*np.sin(theta[i+1,j,tt])*np.sin(phi[i+1,j,tt]) - \
                R*np.sin(theta[i,j,tt])*np.sin(phi[i,j,tt]))**2)

        E = T + U
        energy.append(E)

    theta = theta[1:N+1,1:N+1]
    thdot = thdot[1:N+1,1:N+1]
    phi = phi[1:N+1,1:N+1]
    phdot = phdot[1:N+1,1:N+1]

    print('theta',theta)
    print('thetadot',thdot)
    print('phi',phi)
    print('phidot',phdot)

    print(energy[0],energy[-1])

    return [theta,phi]

#coupled_pendula_3D(10000,5,2)

angles = coupled_pendula_3D(steps,finaltime,N)

def animate_coupled_pend(steps,finaltime,N):

    dt = finaltime/steps
    time = np.linspace(0, finaltime, steps+1)

    x = np.zeros((N,N,steps+1))
    y = np.zeros((N,N,steps+1))
    z = np.zeros((N,N,steps+1))

    for i in range(N):
        for j in range(N):
            x[i,j] = R * np.sin(angles[0][i,j])*np.cos(angles[1][i,j]) + j*2*R
            y[i,j] = R * np.sin(angles[0][i,j])*np.sin(angles[1][i,j]) + i*2*R
            z[i,j] = -R * np.cos(angles[0][i,j])

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal', xlim=(-(2*R+1), (2*N*R+R)),\
        ylim=(-(2*R+1), (2*N*R+R)))
    ax.grid()

    lines = []

    for i in range(N):
        line, = ax.plot([], [], 'o')
        lines.append(line)

    def init():
        for i in range(N):
            lines[i].set_data([], [])

        return lines

    def animate(i):
        thisx = []
        thisy = []
        thisz = []
        for j in range(N):
            for k in range(N):
                thisx.append([])
                thisy.append([])
                thisz.append([])

                thisx[j].append(x[k,j][i])
                thisy[j].append(y[k,j][i])
                thisz[j].append(z[k,j][i])

        for j in range(N):

            lines[j].set_data(thisx[j],thisy[j])

        return lines

    ani = animation.FuncAnimation(fig, animate,interval=1, blit=True, \
                                    init_func=init)

    plt.show()

animate_coupled_pend(steps,finaltime,N)

def threeD(steps,finaltime,N):

    dt = finaltime/steps
    time = np.linspace(0, finaltime, steps+1)

    x = np.zeros((N,N,steps+1))
    y = np.zeros((N,N,steps+1))
    z = np.zeros((N,N,steps+1))

    for i in range(N):
        for j in range(N):
            x[i,j] = R * np.sin(angles[0][i,j])*np.cos(angles[1][i,j]) + j*2*R
            y[i,j] = R * np.sin(angles[0][i,j])*np.sin(angles[1][i,j]) + i*2*R
            z[i,j] = -R * np.cos(angles[0][i,j])

    fig = plt.figure()
    ax = p3.Axes3D(fig)

    ax.set_xlim3d([-2*R-1, 2*N*R+1])
    ax.set_xlabel('X')

    ax.set_ylim3d([-2*R-1, 2*N*R+1])
    ax.set_ylabel('Y')

    ax.set_zlim3d([-R-1, 1])
    ax.set_zlabel('Z')

    ax.set_title('3D NxN Coupled Pendulums')

    lines = []

    for i in range(N):
        line, = ax.plot([], [], [], 'o')
        lines.append(line)

    def init():
        for i in range(N):
            lines[i].set_data([], [])
            lines[i].set_3d_properties([])

        return lines

    def animate(i):
        thisx = []
        thisy = []
        thisz = []
        for j in range(N):
            for k in range(N):
                thisx.append([])
                thisy.append([])
                thisz.append([])

                thisx[j].append(x[k,j][i])
                thisy[j].append(y[k,j][i])
                thisz[j].append(z[k,j][i])

        for j in range(N):
            #set_data() only sets 2 arrays
            lines[j].set_data(thisx[j],thisy[j])
            #need this command to set third dat set
            lines[j].set_3d_properties(thisz[j])

        return lines

    ani = animation.FuncAnimation(fig, animate, interval=1, blit=True, init_func=init)
    plt.show()

threeD(steps,finaltime,N)
