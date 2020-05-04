
library(data.table)
library(deSolve)

rm (list = ls())

sol_ana = function(lambda0, mu0, t)
	return (data.table(lambda = lambda0*exp(-t), mu = -lambda0*exp(-t) + mu0 + lambda0));

lambda0 = c(0.5, sqrt(2), 2, 3.5, 5873.34);
mu0 = c(2, pi, 3.9, 4, 179.58);

sol_ana(lambda0, mu0, 1)

#### With homemade analytical solution
init = fread("init.txt")
init[, .N]
init = init[lambda != 0]
init[, .N]
plot(init$lambda, type = "l", ylim = c(0, 1.35))

end = fread("end.txt")
end[, .N]
end = end[lambda != 0]
end[, .N]
lines(end$lambda)
s = seq(0, 1, length.out = 2000)
lines(exp(-s)*exp(-1), col = "red")
plot(end$lambda, type = "l")
s = seq(0, 1, length.out = end[, .N])
lines((1 - s)*exp(-0.4*3), col = "red")
