% simulate spatio-temporal hawkes process by first generating background
% points (as a collection); then for each background point: generate 
% offsprings and collect these new points;
% the process is iterated until no offspring is generated.

% location of background points are randomly generated, while the location
% of offsprings are generated with respect to the excitation function g

% excitation function g is separable exponential kernels

function [times, locations] = simulate_sthp(v, alpha, beta, grid_min, grid_max, end_time)
alpha_t = alpha(1);
alpha_s = alpha(2);
beta_t = beta(1);
beta_s = beta(2);

% calculate mean number of offspring m 
locfun = @(x,y) exp(-beta_s * sqrt(x.^2+y.^2));
m = integral2(locfun, grid_min, grid_max, grid_min, grid_max);
m = m * alpha_t * alpha_s * (1-exp(-beta_t * end_time)) / beta_t;

%% generating times of background points
T = -log(rand)/v; n = 0;
bg_t = zeros(1, floor(v*end_time));
while T < end_time, n = n+1;
    bg_t(n) = T; T = T - log(rand)/v;
end
bg_t = nonzeros(bg_t).';
[~, bg_count] = size(bg_t);

bg_loc = zeros(2, bg_count);
est_os_count = ceil(m*bg_count);
os_t_container = zeros(1, est_os_count);
os_loc_container = zeros(2, est_os_count);
os_count = 0;


%% looping over background points to generate locations and initial round of offsprings
for point = 1:bg_count
    % generate location
    bg_loc(1, point) = rand*(grid_max-grid_min) + grid_min;
    bg_loc(2, point) = rand*(grid_max-grid_min) + grid_min;

    % simulate number of offsprings
    n_offspring = poissrnd(m);

    % generate offspring times
    temp = exprnd(beta_t, 1, n_offspring);
    os_t = alpha_t * temp / beta_t + bg_t(point);
    
    % generate offspring locations
    temp = exprnd(beta_s, 2, n_offspring);
    signs = fix(rand(2,n_offspring)+0.5);
    signs(~signs)=-1;
    temp = temp .* signs;
    os_loc = alpha_s * temp / beta_s + [bg_loc(1,point); bg_loc(2,point)];
    
    os_t_container(os_count + 1 : os_count + n_offspring) = os_t;
    os_loc_container(:, os_count + 1 :os_count + n_offspring) = os_loc;
    
    os_count = os_count + n_offspring;
end
% disp('offspring points = ');
% disp(os_count);



%% generate offsprings over iterative collections
dummy_indicator = 1;
current_collection = [os_t_container; os_loc_container];
new_t = [];
new_loc = [];
new_count = 0;

while dummy_indicator
    [~, current_count] = size(current_collection);
    collection_t = current_collection(1,:);
    collection_loc = current_collection(2:3,:);
    next_count = 0;
    
    for point = 1:current_count
        % simulate number of offsprings
        n_offspring = poissrnd(m);

        % generate offspring times
        temp = exprnd(beta_t, 1, n_offspring);
        new_os_t = alpha_t * temp / beta_t + collection_t(point);
        % generate offspring locations
        temp = exprnd(beta_s, 2, n_offspring);
        signs = fix(rand(2,n_offspring)+0.5);
        signs(~signs)=-1;
        temp = temp .* signs;
        new_os_loc = alpha_s * temp / beta_s + collection_loc(:, point);
        
        % store in overall container
        new_t(new_count+1 : new_count + n_offspring) = new_os_t;
        new_loc(:, new_count+1: new_count + n_offspring) = new_os_loc;
        new_count = new_count + n_offspring;
    
        % collect for next iterative round
        next_collection(1,next_count+1 : next_count + n_offspring) = new_os_t;
        next_collection(2:3, next_count+1 : next_count + n_offspring) = new_os_loc;
        next_count = next_count + n_offspring;
    end
    
    % escape criteria
    if isempty(next_collection)
        dummy_indicator = 0;
    end
    
    % set up for next round
    current_collection = next_collection;
    next_collection = [];
end

%% collate and thin
% collate all generated points sort wrt time
data = [bg_t, os_t_container, new_t; bg_loc, os_loc_container, new_loc];
[~, order] = sort(data(1,:));
data = data(:,order);


% discard data outside of time and grid
idx = (data(1,:)<end_time).*(data(2,:)<grid_max & data(2,:)>grid_min).*(data(3,:)<grid_max & data(3,:)>grid_min);
data = data.*idx;
times = nonzeros(data(1,:));
x = nonzeros(data(2,:));
y = nonzeros(data(3,:));
locations = [x,y];


