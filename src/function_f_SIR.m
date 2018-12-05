function [X_P, C_value, X_traj, S_traj, P_traj] = function_f_SIR(MCMC_params, X_F, t_n, t_n_plus_1, theta)

t = t_n:MCMC_params.dt:t_n_plus_1;
time_increments = length(t);

%Copy over parameters
dt = MCMC_params.dt;
mu = theta(1);
gamma = theta(2);
R0_avg = theta(3);
alpha = theta(4);
F_noise = theta(5);
beta = R0_avg * (gamma + mu);

%Initial conditions for state variables
S = zeros(time_increments - 1, MCMC_params.J_particles);
S(1,:) = X_F(1,:);
I = zeros(time_increments - 1, MCMC_params.J_particles);
I(1,:) = X_F(2,:);
P = zeros(time_increments - 1, MCMC_params.J_particles);
P(1,:) = X_F(3,:);
R = zeros(time_increments - 1, MCMC_params.J_particles);
R(1,:) = P(1,:) - S(1,:) - I(1,:);
C = zeros(time_increments - 1, MCMC_params.J_particles); %C is cumulative incidence
C(1,:) = 0;

eta_rand = randn(time_increments-1, MCMC_params.J_particles) ./ sqrt(dt); %random variates to add noise

for i = 1:(time_increments-1) %for each dt time increment 
    
    %beta = R0_avg * (gamma + mu);
    seas_now = t(i)/12 - floor(t(i)/12);
    beta_now = beta * (1 + alpha*cos(2*pi*seas_now));
    
    f_term = F_noise .* beta_now .* (S(i,:) ./ P(i,:)) .* I(i,:);
    dN_SI = (beta_now .* (S(i,:)./ P(i,:)) .* I(i,:).* dt) + (f_term .* eta_rand(i,:) .* dt);
    dN_SI(dN_SI < 0) = 0; %make sure this doesn't go negative!
    dN_IR = gamma .* I(i,:) .* dt;
    dN_SD = mu .* S(i,:) .* dt; %mu*S(i) * MCMC_params.dt;
    dN_ID = mu .* I(i,:) .* dt;%mu*I(i) * MCMC_params.dt;
    dN_RD = mu .* R(i,:) .* dt;
    dN_BS(1:MCMC_params.J_particles) = mu .* P(i,:) * dt; %b*P(i) * MCMC_params.dt;

    dS = dN_BS + -dN_SI + -dN_SD;
    dI = dN_SI + -dN_IR + -dN_ID;
    dR = dN_IR + -dN_RD;
    
    S(i+1,:) = S(i,:) + dS;
    I(i+1,:) = I(i,:) + dI;
    R(i+1,:) = R(i,:) + dR;
    
    P(i+1,:) = S(i+1,:) + I(i+1,:) + R(i+1,:); %THIS NEEDS TO BE CHANGED!
    C(i+1,:) = C(i,:) + dN_SI;
    
end

X_traj = I(1:end,:)';
S_traj = S(1:end,:)';
P_traj = P(1:end,:)';
X_P = [S(end, :); I(end, :); P(end, :)];
C_value = C(end,:);
%neg_inds = find(C_value < 0);
%C_value(neg_inds) = 0;