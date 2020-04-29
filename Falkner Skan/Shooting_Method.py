# -*- coding: utf-8 -*-

"""
 Created on 11/19/18
 
 @author: jstickney7
 
 """

import numpy as np
from Falkner_Skan import falk_skan
from RungeKutta import RK5
import matplotlib.pyplot as plt
from matplotlib import rc

def Shoot_secant(beta, guess):

    #Converge value
    chi = 1

    # Calculate Error Value
    error_0 = 1
    error_1 = 1

    # Counter
    i = 0
    while error_0 > 0.001 and error_1 > 0.001:

        eta_0, Sol_0 = falk_skan(beta, guess[i])
        eta_1, Sol_1 = falk_skan(beta, guess[i+1])

        error_0 = abs(Sol_0[-1, 1] - chi)
        error_1 = abs(Sol_1[-1, 1] - chi)

        if error_0 < 0.001:
            print(f"The f''(0) solution is ", guess[i])
            return guess[i], eta_0, Sol_0

        elif error_1 < 0.001:
            print(f"The f''(0) solution is ", guess[i+1])
            return guess[i+1], eta_1, Sol_1
        else:
            #Secant Method for new guess
            guess_new  = guess[i+1] - (Sol_1[-1,1]-1)*(guess[i+1]-guess[i])/(Sol_1[-1,1]-Sol_0[-1,1])
            guess.append(guess_new)
            print(guess_new)

        i += 1

        if i > 1000:
            return "ran out of iterations"


def Shoot_bisect(beta, guess):

    #Converge value
    chi = 1

    # Calculate midpoint error
    error_c = 1

    # Counter
    i = 0

    while abs(error_c) > 0.001:


        # Bisection Method
        c = (guess[1] + guess[0])/2
        print(guess)
        print(c)
        eta_c, Sol_c = falk_skan(beta, c)
        error_c = Sol_c[-1, 1] - chi
        if error_c > 0.001:
            guess[1] = c
        elif error_c < -0.001:
            guess[0] = c
        else:
            print(f"The f''(0) solution is ", c)
            return c, eta_c, Sol_c

        i += 1

        if i > 1000:
            return "ran out of iterations"



if __name__ == "__main__":
    beta = 0.3
    guess = [.5, 1.5]

    ddf, eta, sol = Shoot_bisect(beta, guess)

    esol = zip(eta,sol)

    df = [[eta[i],sol[i,0],sol[i,1],sol[i,2]] for i in range(len(sol)) if sol[i,1] >= 0.999]

    eta_star = df[0][0]-df[0][1]

    theta_star = (ddf - beta*eta_star)/(1+beta)


    print(f"eta_star is ", eta_star)
    print(f"theta_star is ", theta_star)

    fig = plt.figure(1)
    plt.plot(eta, sol[:,0])
    plt.title(r"f vs $\eta$")


    fig = plt.figure(2)
    plt.plot(eta, sol[:,1])
    plt.title(r"f' vs $\eta$")

    fig = plt.figure(3)
    plt.plot(eta, sol[:,2])
    plt.title(f"f'' vs $\eta$")

    plt.show()