# MD
This code is applycable to analyze the collision process (ion-ion/ion-particle/particle-particle) and transportation.

# Requirement
C++ compiler, C++11

# How to use?
**Step1 Compile**

If you use g++ compiler, you can compile the source code by

g++ -o a.out *cpp -std=c++11

**Step2 Set parameters**

in1:  input file 1
con1: calculation condition 1
in2:  input file 2
con2: calculation condition 2
del_coeff:  Coefficient of Limiting sphere (see our JCP paper)
L:    Collision distance
N:    Calculation number
  
**Step3 Run the simulation**

When you run collision simualtion, you have to specify the 7 options. While, in case of diffusion coefficient calculation, 4 options are required as shown following examle.

Collision simulation              : ./a.out in1 con1 in2 con2 del_coeff L N

Diffusion coefficient calculation : ./a.out in1 con1 L N

**Step4 Analyze the output files**
  
