clear
syms E I x_0 y_0 x_1 x_2 x_3 x_4 x_5 x_6 y_1 y_2 y_3 y_4 y_5 y_6 u x_dot_0 x_dot_1 x_dot_2 x_dot_3 x_dot_4 x_dot_5 x_dot_6 y_dot_0 y_dot_1 y_dot_2 y_dot_3 y_dot_4 y_dot_5 y_dot_6 real
syms l positive

% Sextic Bernstein polynomials for 7 control points
b0 = (1-u)^6;
b1 = 6*u*(1-u)^5;
b2 = 15*u^2*(1-u)^4;
b3 = 20*u^3*(1-u)^3;
b4 = 15*u^4*(1-u)^2;
b5 = 6*u^5*(1-u);
b6 = u^6;

P0 = [0 0];
P1 = [x_1 0];
P2 = [x_2 y_2];
P3 = [x_3 y_3];
P4 = [x_4 y_4];
P5 = [x_5 y_5];
P6 = [x_6 y_6];

% Sextic spline curve
pu = b0*P0 + b1*P1 + b2*P2 + b3*P3 + b4*P4 + b5*P5 + b6*P6;

pux = pu(1);
puy = pu(2);

pux_dot = diff(pux, u);
puy_dot = diff(puy, u);
dpu_du = [pux_dot, puy_dot];

pux_2dot = diff(pux_dot, u);
puy_2dot = diff(puy_dot, u); 

%% kappa

%fernet planar curvature
numerator = (pux_dot * puy_2dot - puy_dot * pux_2dot);%abs required?
denominator = (pux_dot^2 + puy_dot^2)^(3/2);
ku =  simplify(numerator / denominator);

ku2 = ku^2;
%%
kugrad = simplify(gradient(ku2, [x_1 x_2 y_2 x_3 y_3 x_4 y_4 x_5 y_5 x_6 y_6]));
%%
%dpu/du
moddpu_du = simplify(sqrt(pux_dot^2 + puy_dot^2));

srtermsimp = simplify(gradient((1- moddpu_du)^2, [x_1 x_2 y_2 x_3 y_3 x_4 y_4 x_5 y_5 x_6 y_6]));

srternosimp = simplify(gradient(((1- moddpu_du)^2)*moddpu_du, [x_1 x_2 y_2 x_3 y_3 x_4 y_4 x_5 y_5 x_6 y_6]));