function [ out ] = toSolve( t, x )
% indiv in cohort = λ = x(1)
% 			  dbh = μ = x(2)

out = zeros(size(x));
% I am 'cheating': I know the analytical value of s*, so I do not estimate it
% s_star = sqrt(1 - 2*d(t, x(2))*t) - 1;
s_star = sqrt(1 - 2*log(exp(d(t, x(2))*t) + exp(-0.5*x(2)^2 - x(2)))) - 1;

out(1) = -d(t, x(2)) * x(1);
out(2) = v(t, x(2), 0);
out(3) = s_star;

end
