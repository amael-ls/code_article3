
#### Aim of prog: Create the initial condition for the C++ prog
## Initial condition from real data
#
## Initial condition from a continuous function (here a lognormal)
#

#### Load packages and clear memory
library(data.table)

rm(list = ls())

#### Create the initial condition N(0, t) = (1 + s)*Exp[-d*t]
# s_inf = 1388;
# nbCohorts = 3000;
# delta_s = s_inf/nbCohorts;
# sMax_initCond = 500;
# sizes = seq(0, sMax_initCond, delta_s)
# treesIC = data.table(density = 1 + sizes, dbh = sizes)

sizes = c(123, 750, 760)
treesIC = data.table(density = 1 + sizes, dbh = sizes)

#### Write the initial condition
write.table(treesIC, file = "./ic_1.txt", append = FALSE, quote = FALSE, sep = " ",
	eol = "\n", na = "NA", dec = ".", row.names = FALSE,
	col.names = TRUE)
