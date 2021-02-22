# MD
This code is applycable to analyze the collision process (ion-ion/ion-particle/particle-particle) and transportation.



# Requirement
C++ compiler, C++11



# How to use?
**Step1 Compile**

  If you use g++ compiler, just comple the codes by

  g++ -o a.out *cpp -std=c++11


**Step2 Set parameters**

  _in1_:  input file 1
  _con1_: calculation condition 1
  _in2_:  input file 2
  _con2_: calculation condition 2
  _del_coeff_:  Coefficient of Limiting sphere (see our JCP paper)
  _L_:    Collision distance
  _N_:    Calculation number
  
  
**Step3 Run the simulation**

  When you run collision simualtion, you have to specify the 7 options. While, in case of diffusion coefficient calculation, 4 options are required as shown following examle.

  Collision simulation              : ./a.out in1 con1 in2 con2 del_coeff L N

  Diffusion coefficient calculation : ./a.out in1 con1 L N


**Step4 Analyze the output files**
  
