function [void] = get_traces(theta_samples, x)

%Make trace plots for all paramters in theta

theta_new = theta_samples(x, :);

[nrows, ncols] = size(theta_new);

for j = 1:nrows
    subplot(nrows, 1, j), plot(theta_new(j,:))
end


end

