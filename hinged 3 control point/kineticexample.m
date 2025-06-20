
clear
close all
syms t u mu x_1 x_2 y_1 y_2 x_dot_0 x_dot_1 x_dot_2 y_dot_0 y_dot_1 y_dot_2 x_ddot_0 x_ddot_1 x_ddot_2 y_ddot_0 y_ddot_1 y_ddot_2 real 
syms l positive
syms x0(t) y0(t) x1(t) x2(t) x3(t) y1(t) y2(t) y3(t)
assume([x0(t), y0(t), x1(t), x2(t), y1(t), y2(t)], 'real')

%berstein polynomials
b0 = (1-u)^3;
b1 = 3*u*(1-u)^2;
b2 = 3*u^2*(1-u);
b3 = u^3;

P0 = [0 0];
P1 = [x1(t) y1(t)];
P2 = [x2(t) y2(t)];
P3 = [x3(t) y3(t)];

q = [0 x1(t) x2(t) x3(t) 0 y1(t) y2(t) y3(t)];%BE CAREFUL ON HOW TO ARRANGE GEN COORD. BE CONSISTENT WITH A!!!
q_dot = diff(q, t)

%spline curve
pu = b0*P0+b1*P1+b2*P2+b3*P3;

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

ds term 
ds = simplify(sqrt(pux_dot^2 + puy_dot^2))
%evaluate ds at rest position
L = l;%length of the beam

%notice that if we evaluate in the rest position u cancells out and ds can be taken out of the integral
ds_cte = simplify(subs(ds, [x1(t) x2(t) y1(t) y2(t)], [L/2 L 0 0]))

A matrix
%matrix terms at each position
% Ax00t = simplify(b0*b0*ds);
% Ax00t = simplify(int(Ax00t, u));
% Ax00 = simplify(subs(Ax00t, u, 1)-subs(Ax00t, u, 0));
% Ax01 = b0*b1*ds
% Ax02 = b0*b2*ds
% Ax11 = b1*b1*ds
% Ax12 = b1*b2*ds
% Ax22 = b2*b2*ds

%matrix terms simplified (ds out)
A00 = int(b0*b0, u, [0 1]);
A01 = int(b0*b1, u, [0 1]);
A02 = int(b0*b2, u, [0 1]);

A11 = int(b1*b1, u, [0 1]);
A12 = int(b1*b2, u, [0 1]);

A22 = int(b2*b2, u, [0 1]);

% Create the upper triangular part of the matrix
upper_tri = [A00, A01, A02;
             0,   A11, A12;
             0,    0,  A22];

% Construct the full symmetric matrix using the upper triangular part
Ax = upper_tri + upper_tri.' - diag(diag(upper_tri))

%Ax = mu*ds_cte* Ax;%mu to be susbtituted by an actual v
mu_val = mu;
Ax = mu_val*ds_cte* Ax;


A = blkdiag(Ax, Ax);

T simplified expression
T_simplified = simplify(1/2 * (q_dot * A * q_dot.'))
symvar(T_simplified)



dT_dq = gradient(T_simplified, [x1(t) x2(t) y1(t) y2(t)]);%second term
dT_dq = simplify(dT_dq)%0 since now ds is cte so there are not q dependant terms in T, only q_dot

dT_dqdot = gradient(T_simplified, [diff(x1(t), t) diff(x2(t), t) diff(y1(t), t) diff(y2(t), t)]);
dT_dqdot = simplify(dT_dqdot)

d_dT_dqdot_dt = diff(dT_dqdot, t);
d_dT_dqdot_dt = simplify(d_dT_dqdot_dt)

% numeric value evaluation to check if computations are equal in both t dependant and t non dependant code
dT_dq_val = vpa(subs(dT_dq, [diff(x1(t), t) diff(x2(t), t) diff(y1(t), t) diff(y2(t), t)], [5 6 7 8]));
dT_dq_val = vpa(subs(dT_dq_val, [x1(t) x2(t) y1(t) y2(t)], [1 2 3 4]))


dT_dqdot_val = vpa(subs(dT_dqdot, [diff(x1(t), t) diff(x2(t), t) diff(y1(t), t) diff(y2(t), t)], [5 6 7 8]));
dT_dqdot_val = vpa(subs(dT_dqdot_val, [x1(t) x2(t) y1(t) y2(t)], [1 2 3 4]))


% subs t dependant for non dependant variables
% x(t) - x
% diff(x(t), t)--> x_dot
% diff(x(t), t, t)--> x_ddot
dT_dq_not = subs(dT_dq, [diff(x1(t), t), diff(x2(t), t), diff(y1(t), t), diff(y2(t), t)], [x_dot_1, x_dot_2, y_dot_1, y_dot_2]);
dT_dq_not = subs(dT_dq_not, [x1(t), x2(t), y1(t), y2(t)], [x_1, x_2, y_1, y_2])
symvar(dT_dq_not)
dT_dq_val_not = vpa(subs(dT_dq_not, [x_1, x_2, y_1, y_2, x_dot_1, x_dot_2, y_dot_1, y_dot_2], [l, l/2, 0, 0, 0, 0, 0, 0]))


dT_dqdot_not = subs(dT_dqdot, [diff(x1(t), t), diff(x2(t), t), diff(y1(t), t), diff(y2(t), t)], [x_dot_1, x_dot_2, y_dot_1, y_dot_2]);
dT_dqdot_not = subs(dT_dqdot_not, [x1(t), x2(t), y1(t), y2(t)], [x_1, x_2, y_1, y_2])
symvar(dT_dqdot_not)
dT_dqdot_val_not = vpa(subs(dT_dqdot_not, [x_1, x_2, y_1, y_2, x_dot_1, x_dot_2, y_dot_1, y_dot_2], [l, l/2, 0, 0, 0, 0, 0, 0]))


d_dT_dqdot_dt_not = subs(d_dT_dqdot_dt, [diff(x1(t), t, t), diff(x2(t), t, t), diff(y1(t), t, t), diff(y2(t), t, t)], [x_ddot_1, x_ddot_2, y_ddot_1, y_ddot_2]);
d_dT_dqdot_dt_not = subs(d_dT_dqdot_dt_not, [diff(x1(t), t), diff(x2(t), t), diff(y1(t), t), diff(y2(t), t)], [x_dot_1, x_dot_2, y_dot_1, y_dot_2])
symvar(d_dT_dqdot_dt_not)
d_dT_dqdot_dt_not = subs(d_dT_dqdot_dt_not, [x1(t), x2(t), y1(t), y2(t)], [x_1, x_2, y_1, y_2]);



