
%%%% In the case without feedback loop
%% Calculate solution for one cohort: N = 1 + π, dbh = π
% Initial condition:
t0 = 0;
t_end = 1;
s0 = pi;
n0 = 1 + s0

[t, y] = ode45(@toSolve, [t0 t_end], [n0, s0]);

%% Analytical solution, cf without_feedbackloop.txt for the c++ output
% Solution (given initial condition is 1 + s)
lambda = @(t, s) (1 + s) .* exp(-d(t, s) * t);
mu_fct = @(t, mu0) sqrt(4*t + mu0^2 + 2*mu0 + 1) - 1;

y(1, :) == [lambda(t0, s0), mu_fct(t0, s0)]

y(end, 1) - lambda(t_end, mu_fct(t0, s0))
y(end, 2) - mu_fct(t_end, s0)
