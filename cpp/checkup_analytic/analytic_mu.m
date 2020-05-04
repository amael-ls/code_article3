
%%%% Aim of prog: Find analytic solution of mu
%% Declare variables
syms mu(t) s_m;
ode = diff(mu, t) == (2 + sqrt(1 - 2*log(exp(d*t) + exp(-s_m^2/2 - s_m))) - 1)*exp(-(sqrt(1 - 2*log(exp(d*t) + exp(-s_m^2/2 - s_m))) - 1))/(1 + mu);

cond = mu(0) == 2;
mu_sol(t) = dsolve(ode, cond);
