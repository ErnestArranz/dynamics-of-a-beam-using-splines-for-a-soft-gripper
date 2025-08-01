clear 
syms E A x_0 y_0 x_1 x_2 x_3 x_4 y_1 y_2 y_3 y_4 u x_dot_0 x_dot_1 x_dot_2 x_dot_3 x_dot_4 y_dot_0 y_dot_1 y_dot_2 y_dot_3 y_dot_4 real

% Quartic Bernstein polynomials
b0 = (1-u)^4;
b1 = 4*u*(1-u)^3;
b2 = 6*u^2*(1-u)^2;
b3 = 4*u^3*(1-u);
b4 = u^4;

P0 = [0 0];
P1 = [x_1 0];
P2 = [x_2 y_2];
P3 = [x_3 y_3];
P4 = [x_4 y_4];

%spline curve
pu = b0*P0+b1*P1+b2*P2+b3*P3+b4*P4;

pux = pu(1);
puy = pu(2);

pux_dot = diff(pux, u);
puy_dot = diff(puy, u);
dpu_du = [pux_dot, puy_dot];

pux_2dot = diff(pux_dot, u);
puy_2dot = diff(puy_dot, u); 

%% 
% 
%% 
% ds

ds = simplify(sqrt(pux_dot^2 + puy_dot^2));

%% 
% gradient (1-ds)^2

grades = gradient((1 - ds)^2, [x_1 x_2 y_2 x_3 y_3 x_4 y_4]);


