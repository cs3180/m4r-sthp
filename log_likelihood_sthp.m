function [log_ll] = log_likelihood_sthp(times, locations, v, alpha, beta, end_time, grid_min, grid_max)


% stationarity tests
% if min(size(times))==1
% 	if  any(sort(times)~=times)
%     		log_ll = -inf;
%     	return
% 	end
% 	if length(unique(times))<length(times)
%     		log_ll = -inf;
%     	return
% 	end
% 	if max(times>end_time)
%     		log_ll = -inf;
%     	return
%     end
%     if max(locations>grid_size)
%         log_ll = -inf;
%         return
%     end

% unpack vectorised params for evaluation
beta_t = beta(1);
beta_s = beta(2);
alpha_t = alpha(1);
alpha_s = alpha(2);

    
% Vectorised version, time-dependent terms and params:
T_s = zeros(size(times));
tp_diff = diff(times);
[total_events,~] = size(times);
for tp = 1:total_events-1
    T_s(tp+1) = (exp(-beta_t*tp_diff(tp))).*(1+T_s(tp));
end

% location-dependent terms and params:
S_s = zeros(size(times));
raw_loc_diff = diff(locations);
sp_diff = sqrt(raw_loc_diff(:,1).^2 + raw_loc_diff(:,2).^2);
for sp = 1:(total_events-1)
    S_s(sp+1) = (exp(-beta_s*sp_diff(sp))).*(1+S_s(sp));
end

% combine and evaluate log-likelihood
f_int = sum(log(v + (alpha_t * T_s) .* (alpha_s * S_s)));

% f_comp = - v*end_time + alpha_s * alpha_t / beta_s / beta_t * sum(exp(-beta_t(end_time-times))-1);

temp_sum = 0;
grid_mesh = linspace(grid_min,grid_max,200);
% calculate grid area delta_g
delta_g = (grid_mesh(2) - grid_mesh(1))^2;

for i = 1:total_events
    int_t = (1 - exp(-beta_t * (end_time - times(i)))) / beta_t;
    % int_t = integral( @(x) exp(-beta_t * (x - times(i))), times(i), end_time);
%     int_s = integral2( @(x, y) exp(-beta_s * ( sqrt((x - locations(1,i)).^2 + (y-locations(2,i).^2)))), grid_min, grid_max, grid_min,grid_max);
    int_s = sum(exp(-beta_s * sqrt((grid_mesh - locations(i,1)).^2 + (grid_mesh-locations(i,2)).^2))*delta_g);
    temp_sum = temp_sum + int_t * int_s;
end

log_ll = f_int - alpha_t * alpha_s * temp_sum - v*end_time*((grid_max-grid_min)^2);
log_ll = -log_ll;
end

