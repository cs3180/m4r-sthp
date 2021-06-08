%% simulate STHPs with the same parameter in bulk

function [data] = simulate_multiple_sthp(v, alpha, beta, grid_min, grid_max, end_time, n_sequences)

% v:                 background parameter mu
% alpha:             excitation rate parameter 
% beta:              decay parameter
% grid_min:          lower bound for spatio axes
% grid_max:          upper bound for spatio axes
%   (spatio limits should be scaled onto the range -1 to 1 for ease of
%   computation. Larger grid would result in long computation time resolving gridmesh used for triple integral of MLE)
% end_time:          maximum time to simulate for (inner-most simulation
% function will only terminate once this is reached.)
% n_sequences:       number of STHP of the same


% initialise container
data = zeros(floor(v*end_time),3,n_sequences);

% loop to generate STHP, which will then occupy the "data" container
for seq = 1:n_sequences
    [times, locations] = simulate_sthp(v, alpha, beta, grid_min, grid_max, end_time);
    new_data = [times,locations];
    [old_length, ~, ~] = size(data);
    [new_length, ~] = size(new_data);
    
    %% update the container wrt length of newly generated sthp to ensure that the resulting data is neat
    if new_length > old_length
        data = cat(1, data, zeros(new_length-old_length,3,n_sequences));
        data(:, :, seq) = new_data;

    
    elseif new_length < old_length
        data_extended = cat(1,new_data, zeros(old_length-new_length,3));
        data(:, :, seq) = data_extended;
 
    
    elseif new_length == old_length
        data(:, :, seq) = new_data;
    end    
    fprintf('finished generating %i th sthp \n', seq);
end

[l, ~, ~] = size(data);
data = permute(data, [1,3,2]);
data = reshape(data, l*n_sequences,3);

% output data to store locally
writematrix(data,'m7b2a1.csv') ;
end


    
    


