function [p_theta, X_traj_sample] = get_Likelihoods(theta, data, dataG, MCMC_params, X_start, X_traj_indexes)

%Set up for SMC
X_traj_new = [];
obsv_times = length(data.y_vals);
X_F = zeros(length(X_start), MCMC_params.J_particles, obsv_times); %X_F is a 3D matrix for the state of each particle through time
X_F(:,:,1) = repmat(X_start, 1, MCMC_params.J_particles);
w_matrix = zeros(MCMC_params.J_particles, obsv_times-1); %w_matrix is for particle weights
k_matrix = zeros(MCMC_params.J_particles, obsv_times-1); %k_matrix holds particle parent indexes

for n = 1:(obsv_times-1) %for each time interval
    index = n + 1;
    t_n = data.t_vals(index - 1); %start of time interval
    t_n_plus_1 = data.t_vals(index); %end of time interval
    start_index = dataG.indices(n); %start index in dataG for time interval 
    end_index = dataG.indices(n+1); %end index in dataG for time interval
    %Run SMC
    [X_F(:,:,index), w_matrix(:,n), X_traj, k_matrix(:,n)] = run_SMC(MCMC_params, X_F(:,:,index-1), theta, t_n, t_n_plus_1, data.y_vals(index), dataG.lineages(start_index:end_index), dataG.coal_events(start_index:end_index), dataG.dt_ref(start_index:end_index), dataG.event_times(start_index:end_index));
    X_traj_new = [X_traj_new, X_traj(:, 2:end)];
end

%Calculate overall likelihood for each particle
log_likelihood = 0;
for i = 1:(length(data.y_vals)-1)
    sum_likelihood_J = 0;
    sum_likelihood_J = sum(w_matrix(:,i)/MCMC_params.J_particles, 1);
    log_likelihood = log_likelihood + log(sum_likelihood_J);
end 
p_theta = log_likelihood; %prob of y_t given theta

%Track one lineage's trajectory to sample from posterior of the latent
%variable X_I
X_traj_1 = repmat(X_start(2), 1, MCMC_params.J_particles);
X_traj_final = [X_traj_1' X_traj_new];
[X_traj_sample] = get_sample_traj(k_matrix, X_traj_final, MCMC_params, X_traj_indexes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate effective sampling size (ESS) to check for particle degeneracy
% W = sum(w_matrix, 1);
% for i = 1:obsv_times-1
%     norm_weights = w_matrix(:, i) ./ W(i);
%     sum_sqr_weights = sum(norm_weights.^2);
%     ESS(i) = 1/sum_sqr_weights;
% end

%View particle trajectories
%I_traj = reshape(X_F(2, :, :), MCMC_params.J_particles, obsv_times);
%keyboard
%I_traj_var = var(I_traj, 1);
%S_traj = reshape(X_F(1, :, :), MCMC_params.J_particles, obsv_times);
