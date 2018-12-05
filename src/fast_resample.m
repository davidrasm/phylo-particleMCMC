function k = fast_resample(w, MCMC_params)

w_norm = w/sum(w);

r = mnrnd(MCMC_params.J_particles, w_norm); %resample by multinomial sampling with replacement

k = zeros(1, MCMC_params.J_particles);
end_index = 0;
for j = 1:MCMC_params.J_particles
    reps = r(j);
    if r(j) > 0
        curr_index = end_index + 1;
        end_index = curr_index + reps - 1;
        k(curr_index:end_index) = j;
    end
end

k_shuffle = randperm(MCMC_params.J_particles);
k = k(k_shuffle);

