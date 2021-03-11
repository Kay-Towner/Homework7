
#By: Kay Towner

import time
import numpy as np
import math

#L1
def func_L1(r, G=None, M=None, m=None, R=None, omega=None):
    """Function that has the formula for lagrange point L1."""
    return ((G*M)/(r**2))-((G*m)/((R-r)**2))- ((omega**2)*r)

def fprime1(r, G, M, m, R, omega):
    """the derivative of func_L1"""
    return ((-2*G*M)/(r**3)) - ((2*G*m)/((R-r)**3)-(omega**2))
#L2
def func_L2(r, G, M, m, R, omega):
    """Function that has the formula for lagrange point L2."""
    return ((G*M)/(r**2))+((G*m)/((R-r)**2))- ((omega**2)*r)

def fprime2(r, G, M, m, R, omega):
    """the derivative of func_L2"""
    return ((-2*G*M)/(r**3)) + ((2*G*m)/((R-r)**3)-(omega**2))


#a.) 
def newtons_method(func=None, deriv=None, firstguess=None,
                   max_iteration=None):
    """Function to run Newton's method on the functions L1 and L2."""
    guess = firstguess
    for i in range(max_iteration):
        print(i, "guess=", guess)
        guess = guess - func(guess)/fprime(guess)
    return guess

#b.)
def secant_method(func=None, deriv=None, firstguess=None, secondguess=None, max_iteration=None):
    """Function to find the Secant Method."""
    x1 = firstguess
    x2 = secondguess
    dumer = (x1 * func(x2)) - (x2 * func(x1)) 
    denom = func(x2) - func(x1)
    for i in range(max_iteration):
        dumer = (x1 * func(x2)) - (x2 * func(x1)) 
        denom = func(x2) - func(x1)
        sec = dumer / denom
        
    return dumer / denom


if __name__ == "__main__":

    G = 6.674*10**(-11) #m^3 kg^-1 s^2  Newton's gravity constant
    M = 1.989*10**30  #kg   Mass of the Sun
    m = 5.972*10**24  #kg  Mass of the Earth
    R = 1.496*10**11  #m  Distance for the sun to Earth
    omega = 1.991*10**(-7) #s^-1  Angular frequency

    #tests for the funcions L1 and L2:
    test = func_L1(G, M, m, R, omega)
    print("This is the result of r for L1:", test, "m")
    testtwo = func_L2(G, M, m, R, omega)
    print("This is the result of r for L2:", testtwo, "m")


    #Answers for parts a. and b. of part 2. :
    print("Newtons Method for L1:", newtons_method(func=func_L1, deriv=fprime1, firstguess=1*10**11,
                   max_iteration=20))
    print("Newtons Method for L2:", newtons_method(func=func_L2, deriv=fprime2, firstguess=1*10**11,
                   max_iteration=20))

    print("Secant Method for L1:", secant_method(func=func_L1, firstguess=1*10**11,
                                secondguess=1*10**11, max_iteration=20))
    
    print("Secant Method for L2:", secant_method(func=func_L2, firstguess=1*10**11,
                                secondguess=1*10**11, max_iteration=20))





    

