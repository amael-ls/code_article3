#### Aim of prog: Check if the solution is oscillatary or not

#### Load library and clear memory
library(data.table)

rm(list = ls())
graphics.off()
options(max.print = 500)

#### Parameters
path = "../cpp/popDyn/Acer_saccharum_fat-tailed/"
min_size = 0
max_size = 10
iter = 2000

#### Load files and extract for (size, time) couple
## List files
ls_files = list.files(path, pattern = "pd_.*.txt")
n = length(ls_files)

ls_results = data.table(filename = ls_files, density = numeric(n))

for (file in ls_files)
{
	temporary = fread(paste0(path, file))
	if (iter %in% temporary[, unique(iteration)])
	{
		temporary = temporary[min_size < dbh & dbh < max_size& iteration == iter]
		if (temporary[, .N] == 0)
		{
			ls_results[filename == file, density := 0]
		} else {
			ls_results[filename == file, density := temporary[dbh == min(dbh), density]]
		}
	} else {
		ls_results[filename == file, density := 0]
	}
}

range(ls_results[, density])
