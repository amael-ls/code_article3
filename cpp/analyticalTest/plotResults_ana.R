
#### Aim of prog: Analytical solution

#### Load packages
library(data.table)
library(tikzDevice)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

options(max.print = 500)

#### Load results c++
## Path
pathCpp = "../"

## Initial state
init = fread(paste0(pathCpp, "init.txt"))

## Last state
end = fread(paste0(pathCpp, "end.txt"))

#### Function analytical solution
## Analytical solution
ana_sol = function(s, t, d)
	return ((1 + s)*exp(-d*t));

## Analytical dbh ()no feedback loop
ana_dbh = function(s0, t)
	return (sqrt(1 + 4*t +2*s0 +s0^2) - 1)

#### Parameters
d = 0.3; # cf C++ code, file Species.c++
t_end = 5

#### Difference analytical versus numerical solution
nonZeroChorts = end[density > 0, .N]

## Absolute errors
end[1:nonZeroChorts, density] - ana_sol(end[1:nonZeroChorts, dbh], t_end, d)
end[1:nonZeroChorts, dbh] - ana_dbh(init[1:nonZeroChorts, dbh], t_end)

## Relative errors
(end[1:nonZeroChorts, density] - ana_sol(end[1:nonZeroChorts, dbh], t_end, d))/abs(ana_sol(end[1:nonZeroChorts, dbh], t_end, d))
(end[1:nonZeroChorts, dbh] - ana_dbh(init[1:nonZeroChorts, dbh], t_end))/abs(ana_dbh(init[1:nonZeroChorts, dbh], t_end))