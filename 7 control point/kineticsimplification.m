clear;
close all;
syms t u mu x_1 x_2 x_3 x_4 x_5 x_6 y_1 y_2 y_3 y_4 y_5 y_6 x_dot_0 x_dot_1 x_dot_2 x_dot_3 x_dot_4 x_dot_5 x_dot_6 y_dot_0 y_dot_1 y_dot_2 y_dot_3 y_dot_4 y_dot_5 y_dot_6 x_ddot_0 x_ddot_1 x_ddot_2 x_ddot_3 x_ddot_4 x_ddot_5 x_ddot_6 y_ddot_0 y_ddot_1 y_ddot_2 y_ddot_3 y_ddot_4 y_ddot_5 y_ddot_6 real; 
syms l positive;
syms x0(t) y0(t) x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) y1(t) y2(t) y3(t) y4(t) y5(t) y6(t);
assume([x0(t), y0(t), x1(t), x2(t), x3(t), x4(t), x5(t), x6(t), y1(t), y2(t), y3(t), y4(t), y5(t), y6(t)], 'real');

%sextic Bernstein polynomials (degree 6)
b0 = (1-u)^6;
b1 = 6*u*(1-u)^5;
b2 = 15*u^2*(1-u)^4;
b3 = 20*u^3*(1-u)^3;
b4 = 15*u^4*(1-u)^2;
b5 = 6*u^5*(1-u);
b6 = u^6;

P0 = [0 0];
P1 = [x1(t) 0];
P2 = [x2(t) y2(t)];
P3 = [x3(t) y3(t)];
P4 = [x4(t) y4(t)];
P5 = [x5(t) y5(t)];
P6 = [x6(t) y6(t)];

q = [0 x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) 0 0 y2(t) y3(t) y4(t) y5(t) y6(t)];%BE CAREFUL ON HOW TO ARRANGE GEN COORD. BE CONSISTENT WITH A!!!
q_dot = diff(q, t);

%spline curve
pu = b0*P0+b1*P1+b2*P2+b3*P3+b4*P4+b5*P5+b6*P6;

%x and y components
pux = pu(1);
puy = pu(2);

pux_subs = subs(pux, [x0(t) y0(t)], [0 0]);
puy_subs = subs(puy, [x0(t) y0(t)], [0 0]);

%first derivative components
pux_dot_time = diff(pux_subs, t);
puy_dot_time = diff(puy_subs, t);

pu_dot_time = [pux_dot_time puy_dot_time];

pux_dot = diff(pux_subs, u);
puy_dot = diff(puy_subs, u);

pux_dot_subs = subs(pux_dot, [x0(t) y0(t)], [0 0]);
puy_dot_subs = subs(puy_dot, [x0(t) y0(t)], [0 0]);

pu_dot = simplify([pux_dot puy_dot]);
%% 
% ds term 

ds = simplify(sqrt(pux_dot^2 + puy_dot^2));
%evaluate ds at rest position
L = l;%length of the beam

%notice that if we evaluate in the rest position u cancells out and ds can be taken out of the integral
ds_cte = simplify(subs(ds, [x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) y1(t) y2(t) y3(t) y4(t) y5(t) y6(t)], [L/6 2*L/6 3*L/6 4*L/6 5*L/6 L 0 0 0 0 0 0]));
%% 
% A matrix

%matrix terms simplified (ds out)
A00 = int(b0*b0, u, [0 1]);
A01 = int(b0*b1, u, [0 1]);
A02 = int(b0*b2, u, [0 1]);
A03 = int(b0*b3, u, [0 1]);
A04 = int(b0*b4, u, [0 1]);
A05 = int(b0*b5, u, [0 1]);
A06 = int(b0*b6, u, [0 1]);

A11 = int(b1*b1, u, [0 1]);
A12 = int(b1*b2, u, [0 1]);
A13 = int(b1*b3, u, [0 1]);
A14 = int(b1*b4, u, [0 1]);
A15 = int(b1*b5, u, [0 1]);
A16 = int(b1*b6, u, [0 1]);

A22 = int(b2*b2, u, [0 1]);
A23 = int(b2*b3, u, [0 1]);
A24 = int(b2*b4, u, [0 1]);
A25 = int(b2*b5, u, [0 1]);
A26 = int(b2*b6, u, [0 1]);

A33 = int(b3*b3, u, [0 1]);
A34 = int(b3*b4, u, [0 1]);
A35 = int(b3*b5, u, [0 1]);
A36 = int(b3*b6, u, [0 1]);

A44 = int(b4*b4, u, [0 1]);
A45 = int(b4*b5, u, [0 1]);
A46 = int(b4*b6, u, [0 1]);

A55 = int(b5*b5, u, [0 1]);
A56 = int(b5*b6, u, [0 1]);

A66 = int(b6*b6, u, [0 1]);

% Create the upper triangular part of the matrix
upper_tri = [A00, A01, A02, A03, A04, A05, A06;
             0,   A11, A12, A13, A14, A15, A16;
             0,    0,  A22, A23, A24, A25, A26;
             0,    0,   0,  A33, A34, A35, A36;
             0,    0,   0,   0,  A44, A45, A46;
             0,    0,   0,   0,   0,  A55, A56;
             0,    0,   0,   0,   0,   0,  A66];

% Construct the full symmetric matrix using the upper triangular part
Ax = upper_tri + upper_tri.' - diag(diag(upper_tri));

mu_val = mu;
Ax = mu_val*ds_cte* Ax;

A = blkdiag(Ax, Ax);
%% 
% T simplified expression

T_simplified = simplify(1/2 * (q_dot * A * q_dot.'));
symvar(T_simplified);
%%

dT_dq = gradient(T_simplified, [x1(t) x2(t) y2(t) x3(t) y3(t) x4(t) y4(t) x5(t) y5(t) x6(t) y6(t)]);%second term
dT_dq = simplify(dT_dq);%0 since now ds is cte so there are not q dependant terms in T, only q_dot

dT_dqdot = gradient(T_simplified, [diff(x1(t), t) diff(x2(t), t) diff(y2(t), t) diff(x3(t), t) diff(y3(t), t) diff(x4(t), t) diff(y4(t), t) diff(x5(t), t) diff(y5(t), t) diff(x6(t), t) diff(y6(t), t)]);
dT_dqdot = simplify(dT_dqdot);

d_dT_dqdot_dt = diff(dT_dqdot, t);
d_dT_dqdot_dt = simplify(d_dT_dqdot_dt);

% subs t dependant for non dependant variables
% x(t) - x
% diff(x(t), t)--> x_dot
% diff(x(t), t, t)--> x_ddot
dT_dq_not = subs(dT_dq, [diff(x1(t), t), diff(x2(t), t), diff(y2(t), t), diff(x3(t), t), diff(y3(t), t), diff(x4(t), t), diff(y4(t), t), diff(x5(t), t), diff(y5(t), t), diff(x6(t), t), diff(y6(t), t)], [x_dot_1, x_dot_2, y_dot_2, x_dot_3, y_dot_3, x_dot_4, y_dot_4, x_dot_5, y_dot_5, x_dot_6, y_dot_6]);
dT_dq_not = subs(dT_dq_not, [x1(t), x2(t), y2(t), x3(t), y3(t), x4(t), y4(t), x5(t), y5(t), x6(t), y6(t)], [x_1, x_2, y_2, x_3, y_3, x_4, y_4, x_5, y_5, x_6, y_6]);
symvar(dT_dq_not);

dT_dqdot_not = subs(dT_dqdot, [diff(x1(t), t), diff(x2(t), t), diff(y2(t), t), diff(x3(t), t), diff(y3(t), t), diff(x4(t), t), diff(y4(t), t), diff(x5(t), t), diff(y5(t), t), diff(x6(t), t), diff(y6(t), t)], [x_dot_1, x_dot_2, y_dot_2, x_dot_3, y_dot_3, x_dot_4, y_dot_4, x_dot_5, y_dot_5, x_dot_6, y_dot_6]);
dT_dqdot_not = subs(dT_dqdot_not, [x1(t), x2(t), y2(t), x3(t), y3(t), x4(t), y4(t), x5(t), y5(t), x6(t), y6(t)], [x_1, x_2, y_2, x_3, y_3, x_4, y_4, x_5, y_5, x_6, y_6]);
symvar(dT_dqdot_not);

d_dT_dqdot_dt_not = subs(d_dT_dqdot_dt, [diff(x1(t), t, t), diff(x2(t), t, t), diff(y2(t), t, t), diff(x3(t), t, t), diff(y3(t), t, t), diff(x4(t), t, t), diff(y4(t), t, t), diff(x5(t), t, t), diff(y5(t), t, t), diff(x6(t), t, t), diff(y6(t), t, t)], [x_ddot_1, x_ddot_2, y_ddot_2, x_ddot_3, y_ddot_3, x_ddot_4, y_ddot_4, x_ddot_5, y_ddot_5, x_ddot_6, y_ddot_6]);
d_dT_dqdot_dt_not = subs(d_dT_dqdot_dt_not, [diff(x1(t), t), diff(x2(t), t), diff(y2(t), t), diff(x3(t), t), diff(y3(t), t), diff(x4(t), t), diff(y4(t), t), diff(x5(t), t), diff(y5(t), t), diff(x6(t), t), diff(y6(t), t)], [x_dot_1, x_dot_2, y_dot_2, x_dot_3, y_dot_3, x_dot_4, y_dot_4, x_dot_5, y_dot_5, x_dot_6, y_dot_6]);
symvar(d_dT_dqdot_dt_not);
d_dT_dqdot_dt_not = subs(d_dT_dqdot_dt_not, [x1(t), x2(t), y2(t), x3(t), y3(t), x4(t), y4(t), x5(t), y5(t), x6(t), y6(t)], [x_1, x_2, y_2, x_3, y_3, x_4, y_4, x_5, y_5, x_6, y_6]);