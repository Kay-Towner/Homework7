#finished by Kay Towner

import time
import numpy as np
import math
import matplotlib.pyplot as plt

def function_to_relax(x, c):
    """Equation to run over the relaxation method."""
    return 1-np.exp(-1*c*x)

def relaxation(start_guess=1, func_to_relax=None,
               func_keyword_args=None, tolerance=1e-6):
    """Function that computes the root via the fixed point
    (relaxation) method. inputs are a starting guess,
    a function to use, any function arguments, and a tolerance
    to exit the function when successive approximations are
    less than this value"""
    guess = start_guess
    guessprime = func_to_relax(guess, c)
    while np.abs(guess - guessprime) > tolerance:
        guess = guessprime
        guessprime = func_to_relax(guessprime, c)  
    return guessprime

if __name__ == "__main__":
    #1. a: 
    c = 2
    print("The solution for part a is:") 
    print("{:.6f}".format(relaxation(func_to_relax=function_to_relax)))

    print("""b) Modify your program to calculate the solution for values of
R0 from 0 to 5 in steps of 0.01 and make a plot of P as a function of R0.
Save this plot as a png file as an output to your program, with properly
labeled axes and titles and such. You should see a clear transition from a
regime in which P = 0 to a regime of nonzero P. This is another example of a
phase transition. In physics this transition is known as the percolation transition;
in epidemiology it is the epidemic threshold.  (7 pts)""")
    #b.
    #loop over r0 values from 0 to 5
    r0_values = np.arange(0, 5.0, 0.01)
    solutions = []
    for r0 in r0_values:
        c = r0
        answer = relaxation(func_to_relax=function_to_relax, func_keywords_args={'c':r0})
        solutions.append(answer)

    #plot of the solutions
    
        plt.plot(solutions, ro_values)
        plt.title("Relaxation Graph")
        plt.savefig('Rexation_Graph.png')
        plt.xlabel('r0')
        plt.ylabel('P')
        plt.show
        

    #save the output data
    output_textfile = 'problem6_1_data.txt'
    np.savetxt(output_textfile,
               np.array(np.vstack((r0_values, solutions))).T,
               delimiter = ', ', header='R_0, Probability of epidemic',
               fmt = ('%.2f', '%.3e'))
