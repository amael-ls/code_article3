
%%%% Case with feedback loop, cf scans 16 september 2019
%% Calculate solution for one cohort: N = 1 + π, dbh = π
% Initial condition:
t0 = 0;
t_end = 5;
s0 = pi;
n0 = 1 + s0
s_star_0 = sqrt(1 - 2*log(1 + exp(-0.5*s0^2 - s0))) - 1;

% Solve numerically (feedback loop must be calculated at the same time)
[t, y] = ode45(@toSolve, [t0 t_end], [n0, s0, s_star_0]);

%% Analytical solution, cf with_feedbackloop.txt for the c++ output
% Solution (given initial condition is 1 + s), mu_fct was provided by mathematica
% These solution works only for a constant deat rate d!
lambda = @(t, s) (1 + s) .* exp(-d(t, s) * t);
mu_fct = @(t, s, mu0) -1  + sqrt(mu0^2 + 2*mu0 + 1 - 14/d(t, s) + 2/d(t, s)*exp(1 - sqrt(1 - 2*d(t, s)*t)) .* ...
	(4 - 2*d(t, s)*t + 3*sqrt(1 - 2*d(t, s)*t)));
% s_star = @(t, s) sqrt(1 - 2*d(t, s)*t) - 1;

y(1, :) - [lambda(t0, s0), mu_fct(t0, s0, s0), 0]

y(end, 1) - lambda(t_end, mu_fct(t0, s0, s0))
y(end, 2) - mu_fct(t_end, s0, s0)

% Check with the last values of end.txt
% cpp = [30.6025, 4.54866] % With s* set to its analytical value
cpp = [30.6025, 5.06138] % With s* calculated
