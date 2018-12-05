function [theta_samples, MCMC_out] = run_MCMC(X_start, data, dataG, MCMC_params, theta_now)

%X_traj is the trajectory of the latent variable I
X_traj_indexes = 1:1/(MCMC_params.dt * 10):((length(data.t_vals) - 1) * (1/MCMC_params.dt)) + 1;
X_traj_length = length(X_traj_indexes);%(length(data.t_vals) - 1) * (1/MCMC_params.dt) * 12 + 1;
X_traj_samples = zeros(MCMC_params.iterations/MCMC_params.log_steps, X_traj_length);

%%%Get initial likelihoods, samples
[p_theta_now, X_traj_now] = get_Likelihoods(theta_now, data, dataG, MCMC_params, X_start, X_traj_indexes);
theta_samples(:,1) = theta_now;
theta_new = theta_now;
p_samples = zeros(1, MCMC_params.iterations); %p_samples holds the marginal likelihoods
p_samples(1) = p_theta_now;
[X_traj_samples(1, :)] = X_traj_now;

%Set up MCMC sampling
log_m = MCMC_params.log_steps:MCMC_params.log_steps:MCMC_params.iterations;
save_m = MCMC_params.save_steps:MCMC_params.save_steps:MCMC_params.iterations;
counter = 2;
proposals = zeros(length(theta_now), MCMC_params.iterations/MCMC_params.log_steps);
accept = zeros(1, 1);

%%%This is the MCMC part%%%
for m = 2:MCMC_params.iterations
    tic
    MCMC_params.m = m;
    MCMC_step = m
    
    %Parameter proposals
    theta_new = mvnrnd(theta_now, MCMC_params.theta_cov);
    while min(theta_new) < 0 || theta_new(6) > 1 %if any parameter is negative OR if rho is greater than one, propose new theta
        %IMPORTANT: mvnrnd is a Matlab function that Octave does not have
        %by default!!
        theta_new = mvnrnd(theta_now, MCMC_params.theta_cov);
    end
    curr_proposal = theta_new;
    
    %Run SMC to get the marginal likelihood for theta_new
    [p_theta_new, X_traj_new] = get_Likelihoods(theta_new, data, dataG, MCMC_params, X_start, X_traj_indexes);
    p_samples(m) = p_theta_new; %these are the marginal likelihoods
        
    %Accept or reject theta proposal
    a = exp(p_theta_new - p_theta_now);
    z = rand;
    if z < a
        X_traj_now = X_traj_new;
        theta_now = theta_new;
        accept = accept + 1;
        p_theta_now = p_theta_new;
    end
    
    if ismember(m, log_m) %log samples
        theta_samples(:,counter) = theta_now;
        proposals(:, counter) = curr_proposal;
        X_traj_samples(counter, :) = X_traj_now;
        counter = counter + 1;
        if ismember(m, save_m) %save samples
            save('MCMC_temp_out', 'X_traj_samples', 'theta_samples')
        end
    end 
    toc
end

MCMC_out.accept_rate = accept/m; %acceptance rate
MCMC_out.p_samples = p_samples;
MCMC_out.proposals = proposals;
MCMC_out.X_samples = X_traj_samples;










    