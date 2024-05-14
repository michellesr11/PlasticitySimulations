function [t_steps, f_steps] = RungeKutta4(model, numVar, CI, t_steps)
% Runge-Kutta 4th-Order Algorithm
f_steps = zeros(length(t_steps), numVar);
f_steps(1, :) = CI;
h = t_steps(2)-t_steps(1);  % Constant time step
for i = 2:length(t_steps)
    k1 = model(t_steps(i-1), f_steps(i-1, :));
    k2 = model(t_steps(i-1)+h/2, f_steps(i-1, :)+k1*h/2);
    k3 = model(t_steps(i-1)+h/2, f_steps(i-1, :)+k2*h/2);
    k4 = model(t_steps(i-1)+h, f_steps(i-1, :)+k3*h);
    f_steps(i, :) = f_steps(i-1, :)+(k1/6+k2/3+k3/3+k4/6)*h;
end
end