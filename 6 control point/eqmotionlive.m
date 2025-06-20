%prepare workspace
clear
close all

                
syms E I A mu x_1 x_2 x_3 x_4 x_5 y_1 y_2 y_3 y_4 y_5 x_dot_0 x_dot_1 x_dot_2 x_dot_3 x_dot_4 x_dot_5 y_dot_0 y_dot_1 y_dot_2 y_dot_3 y_dot_4 y_dot_5 x_ddot_0 x_ddot_1 x_ddot_2 x_ddot_3 x_ddot_4 x_ddot_5 y_ddot_0 y_ddot_1 y_ddot_2 y_ddot_3 y_ddot_4 y_ddot_5 real
syms l F u positive
syms Qben1 Qben2 Qben3 Qben4 Qben5 Qben6 Qben7 Qben8 Qben9 Qstr1 Qstr2 Qstr3 Qstr4 Qstr5 Qstr6 Qstr7 Qstr8 Qstr9 Fy5 real

%load terms\kinetic.mat %we only load kinetic terms because its the only one qddot dependant
% load terms\qben.mat
% load terms\qstr.mat

%% 
% 

%create external force vector
%check all gradient have the same order x_1, x_2, y_2, x_3, y_3, x_4, y_4, x_5, y_5
Fext = [0 0 0 0 0 0 0 0 Fy5]';%force applied to the tip

%check variable dependencies

%subs constants by numeric values
D = 0.02;
E_val = 206e5;
I_val = (pi/64)*D^4;
rho = 1000;
A_val = (pi*D^2)/4; 
mu_val = rho*A_val;
L_val = 1;

%kinetic terms (already derived)
%
k1 = (5*l*mu*x_ddot_1)/99 + (5*l*mu*x_ddot_2)/132 + (5*l*mu*x_ddot_3)/231 + (25*l*mu*x_ddot_4)/2772 + (l*mu*x_ddot_5)/462;
k2 = (5*l*mu*x_ddot_1)/132 + (10*l*mu*x_ddot_2)/231 + (25*l*mu*x_ddot_3)/693 + (5*l*mu*x_ddot_4)/231 + (l*mu*x_ddot_5)/132;
k3 = (10*l*mu*y_ddot_2)/231 + (25*l*mu*y_ddot_3)/693 + (5*l*mu*y_ddot_4)/231 + (l*mu*y_ddot_5)/132;
k4 = (5*l*mu*x_ddot_1)/231 + (25*l*mu*x_ddot_2)/693 + (10*l*mu*x_ddot_3)/231 + (5*l*mu*x_ddot_4)/132 + (2*l*mu*x_ddot_5)/99;
k5 = (25*l*mu*y_ddot_2)/693 + (10*l*mu*y_ddot_3)/231 + (5*l*mu*y_ddot_4)/132 + (2*l*mu*y_ddot_5)/99;
k6 = (25*l*mu*x_ddot_1)/2772 + (5*l*mu*x_ddot_2)/231 + (5*l*mu*x_ddot_3)/132 + (5*l*mu*x_ddot_4)/99 + (l*mu*x_ddot_5)/22;
k7 =(5*l*mu*y_ddot_2)/231 + (5*l*mu*y_ddot_3)/132 + (5*l*mu*y_ddot_4)/99 + (l*mu*y_ddot_5)/22;
k8 =(l*mu*x_ddot_1)/462 + (l*mu*x_ddot_2)/132 + (2*l*mu*x_ddot_3)/99 + (l*mu*x_ddot_4)/22 + (l*mu*x_ddot_5)/11;
k9 =(l*mu*y_ddot_2)/132 + (2*l*mu*y_ddot_3)/99 + (l*mu*y_ddot_4)/22 + (l*mu*y_ddot_5)/11;

kinterm = [k1 k2 k3 k4 k5 k6 k7 k8 k9]';

benterm = [Qben1 Qben2 Qben3 Qben4 Qben5 Qben6 Qben7 Qben8 Qben9]';
strterm = [Qstr1 Qstr2 Qstr3 Qstr4 Qstr5 Qstr6 Qstr7 Qstr8 Qstr9]';
%%

lagsym = kinterm - 0 - benterm -strterm - Fext ==0;

%isolate xyddot coordinates
eq_motion = solve(lagsym, [x_ddot_1 x_ddot_2 y_ddot_2 x_ddot_3 y_ddot_3 x_ddot_4 y_ddot_4 x_ddot_5 y_ddot_5]); %x_1, x_2, y_2, x_3, y_3, x_4, y_4, x_5, y_5

xdd1 = eq_motion.x_ddot_1;
xdd2 = eq_motion.x_ddot_2;
ydd2 = eq_motion.y_ddot_2;
xdd3 = eq_motion.x_ddot_3;
ydd3 = eq_motion.y_ddot_3;
xdd4 = eq_motion.x_ddot_4;
ydd4 = eq_motion.y_ddot_4;
xdd5 = eq_motion.x_ddot_5;
ydd5 = eq_motion.y_ddot_5;

[xdd1 xdd2 ydd2 xdd3 ydd3 xdd4 ydd4 xdd5 ydd5]'

