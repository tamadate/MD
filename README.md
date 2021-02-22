# MD
This code is applycable to analyze the collision process (ion-ion/ion-particle/particle-particle) and transportation.



# Requirement
C++ compiler, C++11



# How to use?
**Step1 Compile**

  g++ -o a.out *cpp -std=c++11

**Step2 Set parameters**

  Collision simulation: Make 4 files and set 3 calculation prperties -> input file 1 (in1), inpute file 2 (in2), calculation condition file 1 (con1), calculation condition file 2 (con2), Coefficient of limiting sphere radius (del_coeff), collision distance (L) and calculation number (N)
  
  Diffusion coefficient calculation: Make 2 files and set 2 calculation prperties -> input file 1 (in1), calculation condition file 1 (con1), collision distance (L) and calculation number (N)
  
  
**Step3 Run the simulation**

  Collision simulation              : ./a.out in1 con1 in2 con2 del_coeff L N

  Diffusion coefficient calculation : ./a.out in1 con1 L N

**Step4 Analyze the output files**
  
#References
DOI: 10.1063/1.5144772 (JCP), DOI: 10.1039/d0cp03989f (PCCP)
