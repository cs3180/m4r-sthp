% MCEM algorithm for (binned) spatio-temporal hawkes processes

function [mean_estimate, params, fval1, grads, hessians] = MCEM_sthp(data, N_monte_carlo, n_times, init_choice, method, seed, end_time, bin_width, tol)

%% initialise parameters and set up for MCEM algorithm

params = zeros(n_times,5); 

if seed~= 0
    rng(seed)
end

params(1,:) = init_choice; % assign random first guess

grads = [];
hessians = [];

% Initialise counters and tolerance 
j=1;
eps=tol+1;

all_weights_unif = [];
all_weights_seq = [];

  while eps>tol && j<n_times %set some tolerance and do while here
        clear T
        fprintf('Current n_times = %i \n', j);
        %EXPECTATION STEP 
        weights_raw = zeros(N_monte_carlo,1);
        for i=1:N_monte_carlo
            if mod(i,10)==0
                fprintf('Current monte carlo loop: %i \n', i);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            current_params = params(j*p-p+1:j*p,:);

            baseline_est = current_params(:,1);
            excitation_est = current_params(:,2:2+p-1);
            decay_est = current_params(:,2+p:end);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Uniform method %
            %[T{i}, log_density] = generate_uniform_times(data,p,bin_width);  
            %weights_raw(i) = -complete_likelihood(T{i},baseline_est,excitation_est,decay_est,end_time,p)-(log_density);
            if strcmp(method,'unif')
                [T(i,:), log_density(i)] = generate_uniform_times(data,p,bin_width);  
                weights_raw(i) = -complete_likelihood_univariate(T(i,:),baseline_est,excitation_est,decay_est,end_time)-(log_density(i));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Sequential method %
            elseif strcmp(method,'seq')
                [T(i,:), ~, log_density(i)] = disc_time_hp_grid(data, 1, [params(j,1),params(j,2),params(j,3)], bin_width, 1); 
                %Compute the Importance Sampling weights 
                weights_raw_seq(i) = exp((complete_likelihood_univariate(T(i,:),params(j,1),params(j,2),params(j,3),end_time)-log_density(i))); 
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%                 disp('Choose a valid method')
%              return
            end
            








