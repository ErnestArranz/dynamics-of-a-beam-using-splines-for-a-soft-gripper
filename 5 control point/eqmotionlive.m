%prepare workspace
clear
close all

                
syms E I A mu x_1 x_2 x_3 x_4 y_1 y_2 y_3 y_4 x_dot_0 x_dot_1 x_dot_2 x_dot_3 x_dot_4 y_dot_0 y_dot_1 y_dot_2 y_dot_3 y_dot_4 x_ddot_0 x_ddot_1 x_ddot_2 x_ddot_3 x_ddot_4 y_ddot_0 y_ddot_1 y_ddot_2 y_ddot_3 y_ddot_4 real
syms l F u positive
syms Qben1 Qben2 Qben3 Qben4 Qben5 Qben6 Qben7 Qstr1 Qstr2 Qstr3 Qstr4 Qstr5 Qstr6 Qstr7 Fy4 real

%load terms\kinetic.mat %we only load kinetic terms because its the only one qddot dependant
% load terms\qben.mat
% load terms\qstr.mat

%% 
% 

%create external force vector
%check all gradient have the same order x_1, x_2, y_2, x_3, y_3, x_4, y_4
Fext = [0 0 0 0 0 0 Fy4]';%force applied to the tip

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
k1 = ;
k2 = ;
k3 = ;
k4 = (8*l*mu*x_ddot_1)/315 + (l*mu*x_ddot_2)/21 + (4*l*mu*x_ddot_3)/63 + (l*mu*x_ddot_4)/18;
k5 = (l*mu*y_ddot_2)/21 + (4*l*mu*y_ddot_3)/63 + (l*mu*y_ddot_4)/18;
k6 = (l*mu*x_ddot_1)/126 + (l*mu*x_ddot_2)/42 + (l*mu*x_ddot_3)/18 + (l*mu*x_ddot_4)/9; 
k7 = (l*mu*y_ddot_2)/42 + (l*mu*y_ddot_3)/18 + (l*mu*y_ddot_4)/9; 

kinterm = [k1 k2 k3 k4 k5 k6 k7]';

benterm = [Qben1 Qben2 Qben3 Qben4 Qben5 Qben6 Qben7]';
strterm = [Qstr1 Qstr2 Qstr3 Qstr4 Qstr5 Qstr6 Qstr7]';
%%

lagsym = kinterm - 0 - benterm -strterm - Fext ==0;

%isolate xyddot coordinates
eq_motion = solve(lagsym, [x_ddot_1 x_ddot_2 y_ddot_2 x_ddot_3 y_ddot_3 x_ddot_4 y_ddot_4]); %x_1, x_2, y_2, x_3, y_3, x_4, y_4

xdd1 = eq_motion.x_ddot_1;
xdd2 = eq_motion.x_ddot_2;
ydd2 = eq_motion.y_ddot_2;
xdd3 = eq_motion.x_ddot_3;
ydd3 = eq_motion.y_ddot_3;
xdd4 = eq_motion.x_ddot_4;
ydd4 = eq_motion.y_ddot_4;

[xdd1 xdd2 ydd2 xdd3 ydd3 xdd4 ydd4]'

