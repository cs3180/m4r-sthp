function [mle_params, max_likelihood] = complete_likelihood_sthp(times, locations, v, alpha, beta, end_time, grid_min, grid_max)

% times           the time stamps 
% Note, in the univariate case this has been vectorised so it is expecting
% a row vector of times, or a matrix where each row is an MC rep
% In the multivariate case, this will take an Nxp matrix (N:number of bins)
% v, a, b         nu, alpha, beta, parameters of exponential kernel
% end time        max simulation time

% grad = -1e10*ones(size([v,a,b]));
% hessian = -1e10*ones(numel([v,a,b]),numel([v,a,b]));

current_params = [v, alpha, beta];
% options1 = optimoptions('fmincon','MaxIter',2e4,'MaxFunEvals',2e4,'Algo','interior-point','TolFun',1e-10,'TolX',1e-10,'GradObj','off','Hessian','off','Display','off','DerivativeCheck','off');
myfun = @(z) log_likelihood_sthp(times, locations, z(1), z(2:3), z(4:5), end_time, grid_min, grid_max);
% [mle_params, max_likelihood] = fmincon(myfun,current_params,[],[],[],[],[],[],[],options1); 
constrain_A = [0, 0, 0, 0, -1; 0, 0, 0, -1, 0; 0, 0, -1, 0, 0; 0, -1, 0, 0, 0];
constrain_b = [0; 0; 0; 0;];
[mle_params, max_likelihood] = fmincon(myfun, current_params, constrain_A, constrain_b); 


% [mle_params, max_likelihood] = fminsearch(myfun, current_params);


