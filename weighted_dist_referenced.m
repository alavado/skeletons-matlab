% Calculates a weighted dist with respect to a reference set
% consisting of bordering voxels, instead of the background.
% input:
%   patches: surface(s) which DT we will calculate.
%   ref_set: bordering voxels in the pss.
%   barrier_set: voxels that don't propagate distance data.
%   weights: [face edge vertex] distance transform weights.
% output:
%   <face, edge, vertex> - distance transform of the patches.
function dt = weighted_dist_referenced(patches, ref_set, barrier_set1, barrier_set2, weights)
    
    dims = size(patches);
    wf = weights(1);
    we = weights(2);
    wv = weights(3);
    
    % All bordering voxels are face-adjacent to the background.
    dt = wf * ref_set;
    
    % We mark barrier set.
    %dt = dt - barrier_set;
    
    % Weights for preceding values through iterations.
    cap_weights = [wv we wv; we wf we; wv we wv];
    row_weights = [we wf we];
    val_weight = wf;
    
    % Block propagation.
    dt(barrier_set1 == 1) = Inf;
    dt(barrier_set2 == 1) = Inf;
    
    % Forward pass (Medial rep, p. 174).
    inf_count = 0;
    for pass = 1:3
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~patches(i, j, k) || ref_set(i, j, k)
                    continue;
                end
                
                % We have to find the minimum distance
                % amongst sorrounding visited voxels.
                cap = reshape(dt(i - 1, j - 1: j + 1, k - 1: k + 1), 3, 3);
                cap(cap == 0) = Inf;
                cap = cap + cap_weights;
                prev_row = reshape(dt(i, j - 1, k - 1: k + 1), 1, 3);
                prev_row(prev_row == 0) = Inf;
                prev_row = prev_row + row_weights;
                prev_val = dt(i, j, k - 1);
                if prev_val == 0
                    prev_val = Inf;
                else
                    prev_val = prev_val + val_weight;
                end
                f_mask = [cap(:)' prev_row prev_val];
                if pass > 1 && dt(i, j, k) > 0
                    f_mask = [f_mask dt(i, j, k)];
                end
                dt(i, j, k) = min(f_mask);
            end
        end
    end
    
    % Backward pass.
    for i = dims(1) - 1: -1 : 2
        for j = dims(2) - 1: -1: 2
            for k = dims(3) - 1: -1: 2
                if ~patches(i, j, k) || ref_set(i, j, k)
                    continue;
                end
                    
                % Same as above, this time in reverse order.
                cap = reshape(dt(i + 1, j - 1: j + 1, k - 1: k + 1), 3, 3);
                cap(cap == 0) = Inf;
                cap = cap + cap_weights;              
                prev_row = reshape(dt(i, j + 1, k - 1: k + 1), 1, 3);
                prev_row(prev_row == 0) = Inf;
                prev_row = prev_row + row_weights;
                prev_val = dt(i, j, k + 1);
                if prev_val == 0
                    prev_val = Inf;
                else
                    prev_val = prev_val + val_weight;
                end
                b_mask = [cap(:)' prev_row prev_val dt(i, j, k)];
                dt(i, j, k) = min(b_mask);
                if(dt(i, j, k) == Inf)
                   inf_count = inf_count+1;
                end
            end
        end
    end
    fprintf('WDT %d not reached in pass %d\n', inf_count, pass);
    inf_count = 0;
    end
    
    % Propagate distance data to extended junction voxels.
    % Forward pass.
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~barrier_set2(i, j, k)
                    continue;
                end
                
                % We have to find the minimum distance
                % amongst sorrounding visited voxels.
                cap = reshape(dt(i - 1, j - 1: j + 1, k - 1: k + 1), 3, 3);
                cap(cap == 0) = Inf;
                cap = cap + cap_weights;
                prev_row = reshape(dt(i, j - 1, k - 1: k + 1), 1, 3);
                prev_row(prev_row == 0) = Inf;
                prev_row = prev_row + row_weights;
                prev_val = dt(i, j, k - 1);
                if prev_val == 0
                    prev_val = Inf;
                else
                    prev_val = prev_val + val_weight;
                end
                f_mask = [cap(:)' prev_row prev_val];
                if pass > 1 && dt(i, j, k) > 0
                    f_mask = [f_mask dt(i, j, k)];
                end
                dt(i, j, k) = min(f_mask);                
            end
        end
    end
    
    % Backward pass.
    for i = dims(1) - 1: -1 : 2
        for j = dims(2) - 1: -1: 2
            for k = dims(3) - 1: -1: 2
                if ~barrier_set2(i, j, k)
                    continue;
                end
                    
                % Same as above, this time in reverse order.
                cap = reshape(dt(i + 1, j - 1: j + 1, k - 1: k + 1), 3, 3);
                cap(cap == 0) = Inf;
                cap = cap + cap_weights;              
                prev_row = reshape(dt(i, j + 1, k - 1: k + 1), 1, 3);
                prev_row(prev_row == 0) = Inf;
                prev_row = prev_row + row_weights;
                prev_val = dt(i, j, k + 1);
                if prev_val == 0
                    prev_val = Inf;
                else
                    prev_val = prev_val + val_weight;
                end
                b_mask = [cap(:)' prev_row prev_val dt(i, j, k)];
                dt(i, j, k) = min(b_mask);
                if(dt(i, j, k) == Inf)
                   inf_count = inf_count+1;
                end
            end
        end
    end
    
    % Voxels that couldn't be reached.
    dt_dasd = dt(dt ~= Inf);
    dt(dt == Inf) = max(dt_dasd(:)) + wv;
end

