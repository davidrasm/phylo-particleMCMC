function [X_traj_sample] = get_sample_traj(k, X_traj_final, MCMC_params, X_traj_indexes)

random_lineage = ceil(rand*MCMC_params.J_particles);
[n_particles, n_times] = size(k);
path = zeros(1, n_times);
path(end) = random_lineage;
for j = n_times:-1:2
    parent_index = path(j);
    path(j-1) = k(parent_index, j-1);
end

sub_n = 1/MCMC_params.dt;
X_traj_sample = X_traj_final(1,1);
for k = 1:length(path)
    start_index = (k-1)*sub_n + 2;
    end_index = k*sub_n + 1;
    X_traj_new = X_traj_final(path(k), start_index:end_index);
    X_traj_sample = [X_traj_sample X_traj_new];
end

X_traj_sample = X_traj_sample(X_traj_indexes);

    
    
    