function [ out ] = toSolve( t, x )
% indiv in cohort = λ = x(1)
% 			  dbh = μ = x(2)

out = zeros(size(x));

out(1) = -d(t, x(2)) * x(1);
out(2) = v(t, x(2));

end
