%% Clear console and variables
% clear all
% clc

disp('Matlab main.m is starting')

%% Parallel pool cluster
% create a local cluster object
pc = parcluster('local');

% explicitly set the JobStorageLocation to the temp directory that was
% created in the sbatch script submission.sh
pc.JobStorageLocation = strcat('/scratch/amael/', getenv('SLURM_JOB_ID'));

slurm_cpus = str2num(getenv('SLURM_CPUS_PER_TASK'));

parpool(pc, str2num(getenv('SLURM_CPUS_PER_TASK')))

%% Load data
% --- Allometries
C0_C1 = readtable('../createMatlabData/C0_C1.csv');
allometries = readtable('../createMatlabData/purves2007_allometries.csv');

C0_C1.Properties.RowNames = C0_C1.parameter;
allometries.Properties.RowNames = allometries.species;

% --- Species names
currentSpecies = '28728-ACE-RUB';

% --- Species-specific integral bounds (dbh corresponding to 45m height)
integral_bounds = readtable('../createMatlabData/dbh_params.csv');
integral_bounds.Properties.RowNames = integral_bounds.species_id;

%% Growth function to compute size at age_max. Autonomous ODE, but 't' is required by ODE45
growth_fct = @(t, dbh, c0, c1, c2, scmu_g, scsd_g, scmu_dbh, scsd_dbh) ...
	exp(scsd_g * (c0 + ...
	c1*(dbh - scmu_dbh)/scsd_dbh + ...
	c2*((dbh - scmu_dbh)/scsd_dbh)^2) + scmu_g);

mortality_fct = @(dbh, c0, c1, c2, scmu_dbh, scsd_dbh) 1 - ...
	exp(-exp(c0 + c1*(dbh - scmu_dbh)/scsd_dbh + c2*((dbh - scmu_dbh)/scsd_dbh)^2));
	
%% Run
% --- Species-specific parameters and data
disp(['species id: ', currentSpecies])

scalingGrowth = readtable('../createMatlabData/growthScaling.csv');
scalingGrowth.Properties.RowNames = scalingGrowth.species_id;
scalingGrowth = scalingGrowth(currentSpecies, {'mu', 'sd'});

dbh_scalingGrowth = readtable('../createMatlabData/growthDbhScaling.csv');
dbh_scalingGrowth.Properties.RowNames = dbh_scalingGrowth.species_id;
dbh_scalingGrowth = dbh_scalingGrowth(currentSpecies, {'mu', 'sd'});

mu_g = scalingGrowth.mu;
sd_g = scalingGrowth.sd;
mu_dbh_g = dbh_scalingGrowth.mu;
sd_dbh_g = dbh_scalingGrowth.sd;

dbh_scalingMortality = readtable('../createMatlabData/mortalityDbhScaling.csv');
dbh_scalingMortality.Properties.RowNames = dbh_scalingMortality.species_id;
dbh_scalingMortality = dbh_scalingMortality(currentSpecies, {'mu', 'sd'});

mu_dbh_m = dbh_scalingMortality.mu;
sd_dbh_m = dbh_scalingMortality.sd;

climate_over_g = readtable(char(strcat('../R0/Matlab_data/', currentSpecies, '/matlabGrowth_above.csv')));
climate_over_g = climate_over_g(:, {'beta0', 'beta1', 'beta2'});
climate_under_g = readtable(char(strcat('../R0/Matlab_data/', currentSpecies, '/matlabGrowth_below.csv')));
climate_under_g = climate_under_g(:, {'beta0', 'beta1', 'beta2'});

climate_over_m = readtable(char(strcat('../R0/Matlab_data/', currentSpecies, '/matlabMortality_above.csv')));
climate_over_m = climate_over_m(:, {'beta0', 'beta1', 'beta2'});
climate_under_m = readtable(char(strcat('../R0/Matlab_data/', currentSpecies, '/matlabMortality_below.csv')));
climate_under_m = climate_under_m(:, {'beta0', 'beta1', 'beta2'});

s_star = integral_bounds(currentSpecies, 'dbh_star10').dbh_star10;
s_inf = integral_bounds(currentSpecies, 'dbh_inf').dbh_inf;
age_max = integral_bounds(currentSpecies, 'age_max').age_max;
local_s = readtable(char(strcat('../R0/results/', currentSpecies, '/local_s_inf.csv')), 'ReadVariableNames', false);
local_s = local_s.Var1;

allometries = allometries(currentSpecies, {'a', 'b', 'T'});

n = height(climate_over_g);
fec = 0.0071;

% --- Create folder to store the results
if ~exist('./results', 'dir')
	mkdir('./results')
end

if ~exist(strcat('./results/', currentSpecies), 'dir')
	mkdir(strcat('./results/', currentSpecies))
end

% --- Define vector of results to save
averageGrowth = zeros(n, 1);
averageMortality = zeros(n, 1);

disp('parfor loop starting')

tic
parfor (j = 1:n)
	fixef_growth_over = climate_over_g(j, :);
	fixef_mortality_over = climate_over_m(j, :);

	fixef_growth_under = climate_under_g(j, :);
	fixef_mortality_under = climate_under_m(j, :);

	local_s_inf = local_s(j);
	current_s_inf = min(local_s_inf, s_inf);

	% If s_inf < s_star, then integrate up to s_inf only. Not that I set t = 0 (its value does not matter)
	if current_s_inf < s_star
		averageGrowth(j) = 1/current_s_inf * integral( @(x) growth_fct(0, x, fixef_growth_under.beta0, fixef_growth_under.beta1, fixef_growth_under.beta2, mu_g, sd_g, mu_dbh_g, sd_dbh_g), 0, current_s_inf, 'ArrayValued', true);

		averageMortality(j) = 1/current_s_inf * integral( @(x) mortality_fct(x, fixef_mortality_under.beta0, fixef_mortality_under.beta1, fixef_mortality_under.beta2, mu_dbh_m, sd_dbh_m), 0, current_s_inf, 'ArrayValued', true);
	
	else
		% Else, integrate from 0 to s_star and from s_star to s_inf. Not that I set t = 0 (its value does not matter)
		averageGrowth(j) = 1/current_s_inf * integral( @(x) growth_fct(0, x, fixef_growth_under.beta0, fixef_growth_under.beta1, fixef_growth_under.beta2, mu_g, sd_g, mu_dbh_g, sd_dbh_g), 0, s_star, 'ArrayValued', true) + 1/(current_s_inf - s_star) * integral( @(x) growth_fct(0, x, fixef_growth_over.beta0, fixef_growth_over.beta1, fixef_growth_over.beta2, mu_g, sd_g, mu_dbh_g, sd_dbh_g), s_star, current_s_inf, 'ArrayValued', true);
		
		averageMortality(j) = 1/s_star * integral( @(x) mortality_fct(x, fixef_mortality_under.beta0, fixef_mortality_under.beta1, fixef_mortality_under.beta2, mu_dbh_m, sd_dbh_m), 0, s_star, 'ArrayValued', true) + 1/(current_s_inf - s_star) * integral( @(x) mortality_fct(x, fixef_mortality_over.beta0, fixef_mortality_over.beta1, fixef_mortality_over.beta2, mu_dbh_m, sd_dbh_m), s_star, current_s_inf, 'ArrayValued', true);
	end
end
toc
csvwrite(char(strcat('./results/', currentSpecies, '/averageGrowth.csv')), averageGrowth)
csvwrite(char(strcat('./results/', currentSpecies, '/averageMortality.csv')), averageMortality)

delete(gcp);
exit;

% https://www.mathworks.com/help/matlab/ref/rowfun.html
% climate.TestAvg = mean(climate{:,4:end}, 2);
