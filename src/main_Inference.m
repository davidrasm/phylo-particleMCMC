function [] = main_Inference()

%Input files
mockData_file = 'SIR_endemic_sim3_mockData'; %File containing mock time series and true parameters
coalTimes_file = 'SIR_endemic_sim3_coalTimes';  %File containing genealogical data (sampling and coalescence times)
output_file = 'SIR_endemic_sim3_MCMCout'; %MCMC output goes here

%Load data files
load(mockData_file) %time series data
load(coalTimes_file) %genealogical data: coal times and sampling times in months since present
load SIR_covmat_041111.mat %covariance matrix for MCMC proposal density

%Get time series data
data.t_vals = t_data - min(t_data); %observation times (months since t0)
data.y_vals = [mockdata]; %case count data (monthly incidence)
data.P = epi_params.N_init; %population size (constant)

%MCMC set up
MCMC_params.J_particles = 200; %number of particles to use in SMC
MCMC_params.iterations = 100; %MCMC steps
MCMC_params.log_steps = 10; %log MCMC samples every x steps
MCMC_params.save_steps = 1000; %save MCMC samples every x steps
MCMC_params.dt = epi_params.dt; %integration time step
MCMC_params.burn_in = 100; %not important

%Get genealogical data
coal_times = coal_times_back; %coal times in months since present (tEnd)!!
sample_times = sample_times_back; %sample times in months since present (tEnd)!!
obsv_times = (max(t_data) - t_data); %observation times in months since present (tEnd) !!
obsv_times = round(obsv_times);
%dataG is the data structure that contains the information used to
%calculate the likelihood of the genealogy over preset time intervals
[G_lineages, G_coal_events, G_lineages_dt, G_indices, G_dt_ref, event_times] = get_G_events(MCMC_params, sample_times, sample_sizes, coal_times, obsv_times);
dataG.coal_events = G_coal_events; %vector that indicates if event is a coalecence events
dataG.lineages = G_lineages; %vector for number of lineages over time
dataG.indices = G_indices; %starting and ending index for each observation interval
dataG.dt_ref = G_dt_ref; %reference to nearest dt integration time
dataG.event_times = event_times; %all events (sampling, coalescence and observation events)

%View lineages over time plot aligned with time series
%subplot(2,1,1), plot(G_lineages_dt)
%subplot(2,1,2), plot(prevalence)
%keyboard

%Initial conditions (can be inferred)
%X_I = [true_params.initConditions(1); true_params.initConditions(2); true_params.initConditions(3); data.P_init]; %just S, E, I and P
X_I = [epi_params.S_init; epi_params.I_init; data.P];

%Initial parameter values
theta_now(1) = epi_params.mu; %mu - host birth/death rate
theta_now(2) = epi_params.gamma; %gamma - rate of recovery
theta_now(3) = epi_params.R0; %R0_avg
theta_now(4) = epi_params.alpha; %alpha - seasonality   
theta_now(5) = epi_params.F_noise; %F_noise - environmental noise
theta_now(6) = epi_params.rho; %rho - reporting rate
theta_now(7) = epi_params.tau; %tau - observation variance

%Covariances for proposal densities
% theta_var(1) = 0;
% theta_var(2) = 0.03;
% theta_var(3) = 0.04;
% theta_var(4) = 0.0003;
% theta_var(5) = 0.0004;
% theta_var(6) = 0.0001;
% theta_var(7) = 0.05;

MCMC_params.theta_cov = cov_mat; %use diag(theta_var) if no starting covariance matrix 

%Run MCMC
[theta_samples, MCMC_out] = run_MCMC(X_I, data, dataG, MCMC_params, theta_now);

save(output_file, 'theta_samples', 'MCMC_out', 'MCMC_params', 'data', 'dataG')

%View marginal posterior density of any parameter
% [beta_summary] = quantile(theta_samples(10, MCMC_params.burn_in:MCMC_params.iterations), [.025 .5 .975]);
% hist(theta_samples(10,MCMC_params.burn_in:MCMC_params.iterations), 20)
% line([beta_summary(1), beta_summary(1)], [ylim], 'LineStyle','--', 'Color', 'r')
% line([beta_summary(2), beta_summary(2)], [ylim], 'LineStyle','-', 'Color', 'r')
% line([beta_summary(3), beta_summary(3)], [ylim], 'LineStyle','--', 'Color', 'r')


%View posterior densities of X_I (infections through time)
% X_I = MCMC_out.X_samples;
% for i = 1:length(X_I(1,:))
%    upper_traj(i) = quantile(X_I(1:301, i), .975);
%    median_traj(i) = quantile(X_I(1:301, i), .5);
%    lower_traj(i) = quantile(X_I(1:301, i), .025);
% end
% plot(median_traj, '-k')
% hold on
% plot(upper_traj, '--r')
% plot(lower_traj, '--r')


