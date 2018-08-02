Using experimental data from uniaxial tests on stainless steel, a statistical inverse problem is solved to determine the stress-strain model. This test is for a linear elastic model, where stress = E*strain, and E is the Young's modulus.

Relevant files:
main.C          Initializes QUESO environment and calls application
compute.C       This file handles the statistical inverse problem (SIP) for estimating the magnitude 'E' of Young's modulus
likelihood.C    The SIP definition requires a user defined likelihood function
e_plots.m       Plots using output data
 