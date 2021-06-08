function [log_ll] = complete_data_log_likelihood_sthp(times, locations, v, alpha, beta, end_time, grid_min, grid_max)

% unpack vectorised params for evaluation
beta_t = beta(1);
beta_s = beta(2);
alpha_t = alpha(1);
alpha_s = alpha(2);


% conditions for stationarity
if any(beta < 0)
    log_ll = inf;
    return
elseif any(alpha<0)
    log_ll = inf;
    return
end
    
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
for sp = 1:total_events-1
    S_s(sp+1) = (exp(-beta_s*sp_diff(sp))).*(1+S_s(sp));
end

% induced CIF at each point
cif = v + (alpha_t * T_s) .* (alpha_s * S_s);

% Calculate triggering probabilities:
g_uij = zeros(total_events, total_events);
Pr_Uij = zeros(total_events, total_events);
Pr_U0 = v ./ cif;

for i = 1:total_events-1
   for j = 1:i
       temp_loc_diff = sqrt((locations(i+1, 1) - locations(j,1))^2 + (locations(i+1,2) - locations(j,2))^2);
       g_uij(i+1,j) = alpha_s * alpha_t * exp(-beta_s * temp_loc_diff) * exp(-beta_t * (times(i+1)-times(j)));
       Pr_Uij(i+1,j) = g_uij(i+1,j) * log(g_uij(i+1,j)) /cif(i+1);
    end
end

% establish grid mesh for estimating triple integral
temp_sum = 0;
grid_mesh = linspace(grid_min,grid_max, 100);
delta_g = (grid_mesh(2) - grid_mesh(1))^2;

% evaluate (estimate) triple integral
for i = 1:total_events
    int_t = (1 - exp(-beta_t * (end_time - times(i)))) / beta_t;
    int_s = sum(exp(-beta_s * sqrt((grid_mesh - locations(i,1)).^2 + (grid_mesh-locations(i,2)).^2))*delta_g);
    temp_sum = temp_sum + int_t * int_s;
end


% evaluate expected complete-data log-likelihood
bg_int = sum(Pr_U0 * log(v));
triggered_int = sum(Pr_Uij, 'all');

log_ll = bg_int + triggered_int - temp_sum - v*end_time*(grid_max - grid_min)^2;
log_ll = -log_ll;
end

