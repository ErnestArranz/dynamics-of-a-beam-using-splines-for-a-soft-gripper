%prepare workspace
clear
close all

                
syms E I A mu x_1 x_2 x_3 x_4 x_5 x_6 y_1 y_2 y_3 y_4 y_5 y_6 x_dot_0 x_dot_1 x_dot_2 x_dot_3 x_dot_4 x_dot_5 x_dot_6 y_dot_0 y_dot_1 y_dot_2 y_dot_3 y_dot_4 y_dot_5 y_dot_6 x_ddot_0 x_ddot_1 x_ddot_2 x_ddot_3 x_ddot_4 x_ddot_5 x_ddot_6 y_ddot_0 y_ddot_1 y_ddot_2 y_ddot_3 y_ddot_4 y_ddot_5 y_ddot_6 real
syms l F u positive
syms Qben1 Qben2 Qben3 Qben4 Qben5 Qben6 Qben7 Qben8 Qben9 Qben10 Qben11 Qstr1 Qstr2 Qstr3 Qstr4 Qstr5 Qstr6 Qstr7 Qstr8 Qstr9 Qstr10 Qstr11 Fy6 real

%load terms\kinetic.mat %we only load kinetic terms because its the only one qddot dependant
% load terms\qben.mat
% load terms\qstr.mat

%% 
% 

%create external force vector
%check all gradient have the same order x_1, x_2, y_2, x_3, y_3, x_4, y_4, x_5, y_5, x_6, y_6
Fext = [0 0 0 0 0 0 0 0 0 0 Fy6]';%force applied to the tip

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
% Note: These are placeholder terms based on the 6-point model
% For a proper implementation, you need to derive the actual 7-point kinetic terms
k1 =       (6*l*mu*x_ddot_1)/143 + (9*l*mu*x_ddot_2)/286 + (8*l*mu*x_ddot_3)/429 + (5*l*mu*x_ddot_4)/572 + (3*l*mu*x_ddot_5)/1001 + (l*mu*x_ddot_6)/1716;

k2 =   (9*l*mu*x_ddot_1)/286 + (5*l*mu*x_ddot_2)/143 + (25*l*mu*x_ddot_3)/858 + (75*l*mu*x_ddot_4)/4004 + (5*l*mu*x_ddot_5)/572 + (l*mu*x_ddot_6)/429;

k3 =                           (5*l*mu*y_ddot_2)/143 + (25*l*mu*y_ddot_3)/858 + (75*l*mu*y_ddot_4)/4004 + (5*l*mu*y_ddot_5)/572 + (l*mu*y_ddot_6)/429;

k4 = (8*l*mu*x_ddot_1)/429 + (25*l*mu*x_ddot_2)/858 + (100*l*mu*x_ddot_3)/3003 + (25*l*mu*x_ddot_4)/858 + (8*l*mu*x_ddot_5)/429 + (l*mu*x_ddot_6)/143;

k5 =                         (25*l*mu*y_ddot_2)/858 + (100*l*mu*y_ddot_3)/3003 + (25*l*mu*y_ddot_4)/858 + (8*l*mu*y_ddot_5)/429 + (l*mu*y_ddot_6)/143;

k6 = (5*l*mu*x_ddot_1)/572 + (75*l*mu*x_ddot_2)/4004 + (25*l*mu*x_ddot_3)/858 + (5*l*mu*x_ddot_4)/143 + (9*l*mu*x_ddot_5)/286 + (5*l*mu*x_ddot_6)/286;

k7 =                         (75*l*mu*y_ddot_2)/4004 + (25*l*mu*y_ddot_3)/858 + (5*l*mu*y_ddot_4)/143 + (9*l*mu*y_ddot_5)/286 + (5*l*mu*y_ddot_6)/286;

k8 =      (3*l*mu*x_ddot_1)/1001 + (5*l*mu*x_ddot_2)/572 + (8*l*mu*x_ddot_3)/429 + (9*l*mu*x_ddot_4)/286 + (6*l*mu*x_ddot_5)/143 + (l*mu*x_ddot_6)/26;

k9 =                               (5*l*mu*y_ddot_2)/572 + (8*l*mu*y_ddot_3)/429 + (9*l*mu*y_ddot_4)/286 + (6*l*mu*y_ddot_5)/143 + (l*mu*y_ddot_6)/26;

k10 =               (l*mu*x_ddot_1)/1716 + (l*mu*x_ddot_2)/429 + (l*mu*x_ddot_3)/143 + (5*l*mu*x_ddot_4)/286 + (l*mu*x_ddot_5)/26 + (l*mu*x_ddot_6)/13;

k11 =                                      (l*mu*y_ddot_2)/429 + (l*mu*y_ddot_3)/143 + (5*l*mu*y_ddot_4)/286 + (l*mu*y_ddot_5)/26 + (l*mu*y_ddot_6)/13;


kinterm = [k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11]';

benterm = [Qben1 Qben2 Qben3 Qben4 Qben5 Qben6 Qben7 Qben8 Qben9 Qben10 Qben11]';
strterm = [Qstr1 Qstr2 Qstr3 Qstr4 Qstr5 Qstr6 Qstr7 Qstr8 Qstr9 Qstr10 Qstr11]';
%%

lagsym = kinterm - 0 - benterm -strterm - Fext ==0;

%isolate xyddot coordinates
eq_motion = solve(lagsym, [x_ddot_1 x_ddot_2 y_ddot_2 x_ddot_3 y_ddot_3 x_ddot_4 y_ddot_4 x_ddot_5 y_ddot_5 x_ddot_6 y_ddot_6]); %x_1, x_2, y_2, x_3, y_3, x_4, y_4, x_5, y_5, x_6, y_6

xdd1 = eq_motion.x_ddot_1;
xdd2 = eq_motion.x_ddot_2;
ydd2 = eq_motion.y_ddot_2;
xdd3 = eq_motion.x_ddot_3;
ydd3 = eq_motion.y_ddot_3;
xdd4 = eq_motion.x_ddot_4;
ydd4 = eq_motion.y_ddot_4;
xdd5 = eq_motion.x_ddot_5;
ydd5 = eq_motion.y_ddot_5;
xdd6 = eq_motion.x_ddot_6;
ydd6 = eq_motion.y_ddot_6;

[xdd1 xdd2 ydd2 xdd3 ydd3 xdd4 ydd4 xdd5 ydd5 xdd6 ydd6]'

