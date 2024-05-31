function [observed_diff, p_value] = stat_test(x, y, num_permutations)
    if nargin < 3
        num_permutations = 1000; % Default number of permutation
    end

    % Compute the observed median difference
    observed_diff = abs(median(y) - median(x));
    
    % Combine the two data groups
    combined = [x; y];
    
    % Counter
    count = 0;
    
    % Permutations
    for i = 1:num_permutations
        % Label permutation
        permuted_indices = randperm(length(combined));
        permuted_x = combined(permuted_indices(1:length(x)));
        permuted_y = combined(permuted_indices(length(x)+1:end));
        
        % Compute the new difference
        new_diff = abs(median(permuted_y) - median(permuted_x));
        
        
        % Compare the new difference with the observed one
        if new_diff >= observed_diff
            count = count + 1;
        end
    end
    
    % Calculer la valeur p
    p_value = count / num_permutations;
end



