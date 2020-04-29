# -*- coding: utf-8 -*-

"""
 Created on 10/25/18
 
 @author: Styxny
 
 """

import numpy as np
import  matplotlib.pyplot as plt
from RungeKutta import RK4, RK5, DormandPrince


def falk_skan(beta,guess):

    #Initial Values/ Tuning Values
    initial = 0.001
    k = 1.001


    # # Eta Array
    def generate_eta(initial, k):
        eta = [initial]
        value = initial
        while value < 6.5:
            value *= k
            eta.append(value)


        return eta
    eta = generate_eta(initial, k)
    # eta = np.arange(0.001,6,0.001)

    #Define Equations
    def falk_skan( Y, eta, beta):
        return np.array([Y[1], Y[2], -1*(Y[0]*Y[2]+beta*(1-Y[1]**2))], dtype=np.float64)


    # Initialize lists to hold derivative values and eta values
    Sol_0 = [0,0,guess]

    Sol, i = DormandPrince(falk_skan,Sol_0,eta,beta)
    # Sol, i = RK5(falk_skan, Sol_0, eta, beta)
    # Sol, i = RK4(falk_skan, Sol_0, eta, beta)

    return eta[:i], Sol[:i,:]


if __name__ == "__main__":

    beta = 1.0
    ff0 = 1.3125
    eta, z = falk_skan(beta, ff0)

    print(z[-1,:])
    fig1= plt.figure(1)
    plt.plot(eta, z[:,1])
    plt.title("f' vs eta")
    # plt.xlim((0,10))
    # plt.ylim((0,1.1))

    fig2= plt.figure(2)
    plt.plot(eta, z[:, 2])
    plt.title("f'' vs eta")

    fig3 = plt.figure(3)
    plt.plot(eta, z[:,0])
    plt.title("f vs eta")
    plt.show()
