
#### Aim of prog: Analytical solution

#### Load packages
library(data.table)
library(tikzDevice)
library(stringi)

#### Clear memory and graphs
rm(list = ls())
graphics.off()

options(max.print = 500)

#### Tool function
read_CppResults = function(filename, nbCohorts)
{
	results = fread(filename)
	liszt = vector(mode = "list", length = 6)
	liszt[1:3] = unlist(results[1:nbCohorts, density])
	liszt[4:6] = unlist(results[1:nbCohorts, dbh])
	return(liszt)
}

#### Load results c++
## Path
pathCpp = "../"

## Initial state
init = fread(paste0(pathCpp, "init.txt"))

## Last state results
nbPoints = c(31, 101, 151, 201, 301, 701, 1201, 1801, 2501, 3001)
n = length(nbPoints)
results = data.table(nbPoints = rep(nbPoints, 2), method = character(2*n),
	density1 = numeric(2*n), density2 = numeric(2*n), density3 = numeric(2*n),
	dbh1 = numeric(2*n), dbh2 = numeric(2*n), dbh3 = numeric(2*n))

#### Function analytical solution
## Analytical dbh ()no feedback loop
ana_dbh = function(s0, t)
	return (sqrt(1 + 4*t +2*s0 +s0^2) - 1)

## Analytical solution
ana_sol = function(s0, t, d)
	return ((1 + ana_dbh(s0, t))*exp(-d*t));

#### Parameters
d = 0.3; # cf C++ code, file Species.c++
t_end = 5
nonZeroChorts = 3

results[, delta_t := t_end/(nbPoints - 1)]

#### Run
## Current folder
R_folder = getwd()
C_folder = stri_sub(R_folder, to = stri_locate_last(R_folder, regex = "/")[1] - 1)
subfolder = stri_sub(R_folder, from = stri_locate_last(R_folder, regex = paste0(C_folder, "/"))[2] + 1)

setwd(C_folder)

for (i in 1:n)
{
	system(paste0("sed -i.tmp '2s/.*/n_t: ", results[i, nbPoints], "/' ", subfolder, "/simulationParameters.txt"))
	system(paste0("./euler ", subfolder, "/simulationParameters.txt"))
	results[i, method := "euler"]
	results[i, 3:8 := read_CppResults("end.txt", nonZeroChorts)]
	system(paste0("./rk4 ", subfolder, "/simulationParameters.txt"))
	results[i + n, method := "rk4"]
	results[i + n, 3:8 := read_CppResults("end.txt", nonZeroChorts)]
}

#### Difference analytical versus numerical solution
ls_names = c(paste0("density", 1:nonZeroChorts), paste0("dbh", 1:nonZeroChorts))
diff = data.table(density1 = numeric(2*n), density2 = numeric(2*n), density3 = numeric(2*n),
	dbh1 = numeric(2*n), dbh2 = numeric(2*n), dbh3 = numeric(2*n))

diff[, density1 := abs(results[, density1] - ana_sol(init[1, dbh], t_end, d))]
diff[, density2 := abs(results[, density2] - ana_sol(init[2, dbh], t_end, d))]
diff[, density3 := abs(results[, density3] - ana_sol(init[3, dbh], t_end, d))]

diff[, dbh1 := abs(results[, dbh1] - ana_dbh(init[1, dbh], t_end))]
diff[, dbh2 := abs(results[, dbh2] - ana_dbh(init[2, dbh], t_end))]
diff[, dbh3 := abs(results[, dbh3] - ana_dbh(init[3, dbh], t_end))]

diff[, ana1 := abs((results[, density1] - ana_sol(init[1, dbh], t_end, d))/ana_sol(init[1, dbh], t_end, d))]
diff[, ana2 := abs((results[, density2] - ana_sol(init[2, dbh], t_end, d))/ana_sol(init[2, dbh], t_end, d))]
diff[, ana3 := abs((results[, density3] - ana_sol(init[3, dbh], t_end, d))/ana_sol(init[3, dbh], t_end, d))]

#### Plots
plot(log10(results[1:n, delta_t]), log10(diff[1:n, density1]), type = "l", col = "#3A28C8")
points(log10(results[1:n, delta_t]), log10(diff[1:n, density1]), pch = 20, col = "#3A28C8")
lines(log10(results[1:n + n, delta_t]), log10(diff[1:n + n, density1]), col = "#9B3D22")
points(log10(results[1:n + n, delta_t]), log10(diff[1:n + n, density1]), pch = 20, col = "#9B3D22")

plot(log10(results[1:n, delta_t]), log10(diff[1:n, ana1]), type = "l", col = "#3A28C8")
points(log10(results[1:n, delta_t]), log10(diff[1:n, ana1]), pch = 20, col = "#3A28C8")
lines(log10(results[1:n + n, delta_t]), log10(diff[1:n + n, ana1]), col = "#9B3D22")
points(log10(results[1:n + n, delta_t]), log10(diff[1:n + n, ana1]), pch = 20, col = "#9B3D22")

slopes_euler = (diff[1:(n-1), log10(ana1)] - diff[2:n, log10(ana1)])/(results[1:(n-1), log10(delta_t)] - results[2:n, log10(delta_t)])
slopes_rk = (diff[2:n + n, log10(ana1)] - diff[1:(n-1) + n, log10(ana1)])/(results[2:n + n, log10(delta_t)] - results[1:(n-1) + n, log10(delta_t)])

# for (i in 1:3)
# {
# 	diff[, ls_names[i] := results[, ls_names[i]]]
# 	diff[, ls_names[i + nonZeroChorts] := results[, ..ls_names[i + nonZeroChorts]]]
# }
