function [data] = simulate_multiple_sthp(v, alpha, beta, grid_min, grid_max, end_time, n_sequences)

data = zeros(floor(v*end_time),3,n_sequences);

for seq = 1:n_sequences
    [times, locations] = simulate_sthp(v, alpha, beta, grid_min, grid_max, end_time);
    new_data = [times,locations];
    [old_length, ~, ~] = size(data);
    [new_length, ~] = size(new_data);
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

writematrix(data,'m7b2a1.csv') ;
end


    
    


