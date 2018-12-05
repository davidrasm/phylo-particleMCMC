function [X_F_next, w_vector, X_traj, k] = run_SMC(MCMC_params, X_F, theta_new, t_n, t_n_plus_1, y_n, lineages, coal_events, dt_ref, event_times)

%Simulate particle trajectories

[X_P, C_value, X_traj, S_traj, P_traj] = function_f_SIR(MCMC_params, X_F, t_n, t_n_plus_1, theta_new);

%Calculate particle weights based on case reports
rho = theta_new(6); %reporting rate
tau = theta_new(7); %obsv. variance
%w_ts = pdf('Normal', y_n, rho * C_value, sqrt(rho * tau * C_value)); %weight by the prob of observing time series data (y_n)
w_ts = normpdf(y_n, rho * C_value, sqrt(rho * tau * C_value)); %this works in matlab and octave, pdf('Normal',..) has no support in Octave
w_ts = w_ts';

%Calculate particle weights based on genealogy
w_G = pdTreeMatrix(MCMC_params, theta_new, lineages, coal_events, event_times, dt_ref, X_traj, S_traj, P_traj);

w = w_ts .* w_G; %w is the joint probability of the time series and genealogy
%w = w_G %just weight based on the genealogy to compute p(G|x1:T, theta)
%w = ts %just weight based on time series to compute p(y1:T|x1:T, theta)

%Deal with particles with zero probability
zeroLocs = find(w == 0);
w(zeroLocs) = Inf;
w(zeroLocs) = min(min(w)); %set all particles with zero probs to smallest nonzero prob
if min(min(w)) == Inf
    w(zeroLocs) = realmin; %set all particles to lowest floating point number in matlab if all particles have zero probability
end

%Resample particles for next time step
w_vector = w;
k = fast_resample(w, MCMC_params); %fast_resample uses multinomial resampling
k = k';
X_F_next = X_P(:,k); %copy over state variables for next time interval
