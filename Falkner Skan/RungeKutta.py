# -*- coding: utf-8 -*-

"""
 Created on 11/19/18
 
 @author: Styxny
 
 """

import numpy as np



def RK4(func, y_1, eta, beta):
    """

    :param func: Vector of the functions
    :param y_1: Initial values for the function and it's derivatives
    :param eta: Vector of independent variable steps
    :return: List of lists where each list gives the function/derivatives at the corresponding independent variable step
    """

    # Initialize returned matrix
    Sol = np.zeros((np.size(eta),np.size(y_1)))

    # Set initial Values
    Sol[0,:] = y_1

    for i, n in enumerate(eta[0:-1]):
        # print(i)
        #Step length
        h = eta[i+1]-eta[i]

        h2 = h/2

        #Predictor Calculation
        k1 = np.asarray(func(Sol[i,:], n, beta))

        k2 = np.asarray(func(Sol[i,:]+ k1*h2, n + h2, beta))

        k3 = np.asarray(func(Sol[i,:]+ k2*h2, n + h2, beta))

        k4 = np.asarray(func(Sol[i,:] + k3, n + h, beta))

        Sol[i+1,:] = Sol[i,:] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

        if Sol[i+1,1] < 0:
            break

        if Sol[i+1,1] > 15:
            break

        if np.abs(Sol[i+1,1] - Sol[i,1]) < 0.001:
            break

    # print(Sol[0, :])
    return Sol, (i + 2)



def RK5(func, y_1, eta, beta):
    """

    :param func: Vector of the functions
    :param y_1: Initial values for the function and it's derivatives
    :param eta: Vector of independent variable steps
    :return: List of lists where each list gives the function/derivatives at the corresponding independent variable step
    """

    # Initialize returned matrix
    Sol = np.zeros((np.size(eta),np.size(y_1)))

    # Set initial Values
    Sol[0,:] = y_1

    for i, n in enumerate(eta[0:-1]):
        # print(i)
        #Step length
        h = eta[i+1]-eta[i]

        #Predictor Calculation
        k1 = np.asarray(func(Sol[i,:], n, beta))

        k2 = np.asarray(func(Sol[i,:]+ k1*(h/4), n + (h/4), beta))

        k3 = np.asarray(func(Sol[i,:]+ (h/8)*k1 + (h/8)*k2, n + (h/4), beta))

        k4 = np.asarray(func(Sol[i,:] - 0.5*h*k2 + h*k3, n + h/2, beta))

        k5 = np.asarray(func(Sol[i,:] + (3/16) *h* k1 + (9/16) * h *k4, n + 0.75*h, beta))

        k6 = np.asarray(func(Sol[i,:] - (3/7)*h*k1 + (2/7)*h*k2 + (12/7)*h*k3 - (12/7)*h*k4 + (8/7)*h*k5, n + h, beta))

        Sol[i+1,:] = Sol[i,:] + (h/90) * (7*k1 + 32*k3 + 12*k4 + 32*k5+7*k6)



        if Sol[i+1,1] < 0:
            break

        if Sol[i+1,1] > 2:
            break

        if np.abs(Sol[i+1,1] - Sol[i,1]) < 0.001:
            break

    # print(Sol[0, :])
    return Sol, (i + 2)


def DormandPrince(func, y_1, eta, beta):
    """

     :param func: Vector of the functions
     :param y_1: Initial values for the function and it's derivatives
     :param eta: Vector of independent variable steps
     :return: List of lists where each list gives the function/derivatives at the corresponding independent variable step
     """

    # Initialize returned matrix
    Sol = np.zeros((np.size(eta), np.size(y_1)))

    # Set initial Values
    Sol[0, :] = y_1

    for i, n in enumerate(eta[0:-1]):
        # print(i)
        # Step length
        h = eta[i + 1] - eta[i]

        # Predictor Calculation
        k1 = np.asarray(func(Sol[i, :], n, beta))

        k2 = np.asarray(func(Sol[i, :] + k1 * (h / 5), n + (h / 5), beta))

        k3 = np.asarray(func(Sol[i, :] + (h*3 / 40) * k1 + (h*9 / 40) * k2, n + h*0.3, beta))

        k4 = np.asarray(func(Sol[i, :] + (44*h/45)*k1 - (56/15) * h * k2 + (32/9)*h * k3, n + h*.8, beta))

        k5 = np.asarray(func(Sol[i, :] + (19372/6561) * h * k1 - (25360/2187) * h * k2 + (64448/6561)*h*k3 - (212/729)*h*k4, n + (8/9) * h, beta))

        k6 = np.asarray(func(
            Sol[i, :] + (9017/3168) * h * k1 - (355 / 33) * h * k2 + (46732 / 5247) * h * k3 + (49 / 176) * h * k4 - (5103 / 18656) * h * k5,
            n + h, beta))



        Sol[i + 1, :] = Sol[i, :] + (h) * ((35/384) * k1 + (500/1113) * k3 + (125/192) * k4 - (2187/6784) * k5 + (11/84) * k6)

        # if np.abs(Sol[i + 1:1] - Sol[i:1]) < 0.001:
        #     break


        if Sol[i+1,1] < 0:
            break

        if Sol[i+1,1] > 1.01:
            break

        # if np.abs(Sol[i+1,1] - Sol[i,1]) < 0.001:
        #     break


    # print(Sol[0, :])
    return Sol, (i + 2)
