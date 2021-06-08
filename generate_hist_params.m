format long g
v = 7.01;
alpha = [1.12, 0.981];
beta = [2.012, 1.989];
end_time = 50;
grid_min = -1;
grid_max = 1;
n_sequences = 200;

params_container = zeros(n_sequences,5);
likelihood_container = zeros(n_sequences);

[total_data, ~] = size(a);
seq_length = total_data / n_sequences;

for seq = 1:n_sequences
    data = a((seq-1)*seq_length+1 : seq*seq_length,:);
    times = nonzeros(data(:, 1));
    locations = nonzeros(data(:, 2:3));
    temp = size(locations);
    locations = reshape(locations, temp(1)/2, 2);
    
    % fmincon on loglikelihood
    [mle_params, max_likelihood] = mle_complete_data(times, locations, v, alpha, beta, end_time, grid_min, grid_max);
    params_container(seq,:) = mle_params;
    likelihood_container(seq) = max_likelihood;
    
    fprintf('finished operation for %i th sthp \n', seq);
end

histogram(params_container(:,1));

    
    