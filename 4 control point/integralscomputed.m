syms E I x_0 y_0 x_1 x_2 x_3 y_1 y_2 y_3 u  x_dot_0 x_dot_1 x_dot_2 x_dot_3 y_dot_0 y_dot_1 y_dot_2 y_dot_3 real
syms l positive
L=1;

subs(kugrad(1), [x_1 x_2 y_2 x_3 y_3], [L/3, 2*L/3, 0, L, 0])
