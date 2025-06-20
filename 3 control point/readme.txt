Last simulations obtained.
 
stretching contains the computations for the generalized elastic forces due to stretching.
petentialdssimplif same but for bending
Kinetic contains the kinetic energy term of the lagrange equations

eqmotion live obtains the equations of motion. assembles the lagrange system and isolates qddot variables. 

Mind that all terms should be derived/assembled consistently (x1 x2 y1 y2) in this order.

dynamics contains all the dinamics solved. all integrals are computed numerically due to singularities in the analytical expressions. (by evaluating the particular generalized coordinates in every time step singularities dont appear). 

splineanimation animates the spline behaviour after running dynamics.