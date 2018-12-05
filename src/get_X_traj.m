function [X_traj_sample] = get_X_traj(theta, data, dataG, MCMC_params, X_start, X_traj_indexes)

X_traj_new = [];
obsv_times = length(data.y_vals);
X_F = zeros(length(X_start), MCMC_params.J_particles, obsv_times); %state variables for forward simulations
X_F(:,:,1) = repmat(X_start, 1, MCMC_params.J_particles);
w_matrix = zeros(MCMC_params.J_particles, obsv_times-1);
k_matrix = zeros(MCMC_params.J_particles, obsv_times-1);
%Monthly_Incidence = zeros(MCMC_params.J_particles, obsv_times - 1);

for n = 1:(obsv_times-1)
    index = n + 1;
    [X_F(:,:,index), w_matrix(:,n), X_traj, C_value, k_matrix(:,n)] = run_SMC_X_traj(MCMC_params, X_F(:,:,index-1), theta, data.t_vals(index - 1), data.t_vals(index), data.y_vals(index), data.births(index), data.deaths(index), dataG.G_events{index-1}, dataG.lineages{index-1}, dataG.coal_ind{index-1}, dataG.subtimes{index-1}, dataG.delta_gen);
    X_traj_new = [X_traj_new, X_traj(:, 2:end)];
    %Monthly_Incidence(:, n) = C_value;
end

% W = sum(w_matrix, 1);
% for i = 1:obsv_times-1
%     norm_weights = w_matrix(:, i) ./ W(i);
%     sum_sqr_weights = sum(norm_weights.^2);
%     ESS(i) = 1/sum_sqr_weights;
% end
% keyboard

%%%%%%%Calculate the total likelihood%%%%%%
%log_likelihood = 0;
%for i = 1:(length(data.y_vals)-1)
%    sum_likelihood_J = 0;
%    sum_likelihood_J = sum(w_matrix(:,i)/MCMC_params.J_particles, 1);
%    log_likelihood = log_likelihood + log(sum_likelihood_J);
%end 
%p_theta = log_likelihood; %prob of y_t given theta

%%%Track one lineage's trajectory to sample from X_I%%%
X_traj_1 = repmat(X_start(3), 1, MCMC_params.J_particles);
X_traj_final = [X_traj_1' X_traj_new];
X_traj_sample = get_sample_traj(k_matrix, X_traj_final, MCMC_params, X_traj_indexes);


%%%%View and plot particle trajectories%%%%
%I_traj = reshape(X_F(3, :, :), MCMC_params.J_particles, obsv_times);
%I_traj_var = var(I_traj, 1);
%S_traj = reshape(X_F(1, :, :), MCMC_params.J_particles, obsv_times);
%E_traj = reshape(X_F(2, :, :), MCMC_params.J_particles, obsv_times);
%for i = 1:(obsv_times-1)
%    upper_traj(i) = quantile(Monthly_Incidence(:, i), .99);
%    median_traj(i) = quantile(Monthly_Incidence(:, i), .5);
%    lower_traj(i) = quantile(Monthly_Incidence(:, i), .01);
%end
%subplot(1,2,1),plot(Monthly_Incidence(:,1:120)')
%subplot(1,2,2),hold on,plot(upper_traj(1:120), '--r'),plot(lower_traj(1:120), '--r'),plot(data.y_vals(2:121), '-k', 'LineWidth', 2)
%subplot(1,2,2),plot(lower_traj(1:120), '--r')
%subplot(2,1,2),plot(median_traj(1:120), '-k')
%subplot(1,2,2),plot(data.y_vals(2:121), '-k', 'LineWidth', 2)