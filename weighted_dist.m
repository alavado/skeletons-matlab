% V's <wf, we, wv> weighted distance transform.
% Ref: G. Borgefors 1996
function dt = weighted_dist(V, weights)
    
    pad = 1;
    V = padarray(V, [pad pad pad]);
    dims = size(V);
    wf = weights(1);
    we = weights(2);
    wv = weights(3);
    dt = V;
    
    cap_weights = [wv we wv we wf we wv we wv];
    row_weights = [we wf we];
    val_weight = wf;
    
    % Forward pass.
    for i = 2: dims(1) - 1
        for j = 2: dims(2) - 1
            for k = 2: dims(3) - 1
                if ~V(i, j, k)
                    continue;
                end
                cap = dt(i - 1, j - 1: j + 1, k - 1: k + 1);
                cap = cap(:)' + cap_weights;
                prev_row = dt(i, j - 1, k - 1: k + 1);
                prev_row = prev_row(:)' + row_weights;
                prev_val = dt(i, j, k - 1) + val_weight;
                dt(i, j, k) = min([cap prev_row prev_val]);
            end
        end
    end
    
    % Backward pass.
    for i = dims(1) - 1: -1 : 2
        for j = dims(2) - 1: -1: 2
            for k = dims(3) - 1: -1: 2
                if ~V(i, j, k)
                    continue;
                end
                
                % Same as above, this time in reverse order.
                cap = dt(i + 1, j - 1: j + 1, k - 1: k + 1);
                cap = cap(:)' + cap_weights;                    
                prev_row = dt(i, j + 1, k - 1: k + 1);
                prev_row = prev_row(:)' + row_weights;
                prev_val = dt(i, j, k + 1) + val_weight;
                dt(i, j, k) = min([cap prev_row prev_val dt(i, j, k)]);
            end
        end
    end
    
%     % Weights for preceding values through iterations.
%     cap_weights = [wv we wv; we wf we; wv we wv];
%     row_weights = [we wf we];
%     val_weight = wf;
%     
%     % Forward pass (Medial rep, p. 174).
%     for i = 2: dims(1) - 1
%         for j = 2: dims(2) - 1
%             for k = 2: dims(3) - 1
%                 if ~V(i, j, k)
%                     continue;
%                 end
%                     
%                 % We have to find the minimum distance
%                 % amongst sorrounding visited voxels.
%                 cap = reshape(dt(i - 1, j - 1: j + 1, k - 1: k + 1), 3, 3);
%                 cap = cap + cap_weights;
%                 prev_row = reshape(dt(i, j - 1, k - 1: k + 1), 1, 3);
%                 prev_row = prev_row + row_weights;
%                 prev_val = dt(i, j, k - 1) + val_weight;
%                 f_mask = [cap(:)' prev_row prev_val];
%                 dt(i, j, k) = min(f_mask);
%             end
%         end
%     end
%     
%     % Backward pass.
%     for i = dims(1) - 1: -1 : 2
%         for j = dims(2) - 1: -1: 2
%             for k = dims(3) - 1: -1: 2
%                 if ~V(i, j, k)
%                     continue;
%                 end
%                     
%                 % Same as above, this time in reverse order.
%                 cap = reshape(dt(i + 1, j - 1: j + 1, k - 1: k + 1), 3, 3);
%                 cap = cap + cap_weights;                    
%                 prev_row = reshape(dt(i, j + 1, k - 1: k + 1), 1, 3);
%                 prev_row = prev_row + row_weights;
%                 prev_val = dt(i, j, k + 1) + val_weight;
%                 b_mask = [cap(:)' prev_row prev_val dt(i, j, k)];
%                 dt(i, j, k) = min(b_mask);
%             end
%         end
%     end
%     
    % By definition, distance is infinite outside the shape.
   % dt(dt == 0) = Inf;
   % maybe not

    % We perform label substitution 3 -> 1 for T = 4,
    % according to Arcelli and Sanniti di Baja '88
    dt(dt == wf) = 1;
   
    % Unpad.
    dt = dt(1 + pad: end - pad, 1 + pad: end - pad, 1 + pad: end - pad);
end

