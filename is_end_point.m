% End point test for 2D and 3D binary images.
% input: N-dimensional binary matrix, cell array of integers (point coords.).
% output: 1 if point is end point in V, 0 otherwise.
function endness = is_end_point(V, point)

    % In 2D, P is an end point if it has a single
    % neighbor or two which are 4-adjacent to one another. 
    % (Siddiqi et al, 2002)
    if numel(point) == 2
        
        % Point coordinates.
        i = point{1};
        j = point{2};
        
        % Must be an object point.
        if ~V(i, j)
            endness = false;
            return;
        end
                   
        % We put all neighbors in a list.
        P_neighbors = [];
        for n_i = -1: 1: 1
            for n_j = -1: 1: 1
                if (n_i ~= 0 || n_j ~= 0) && V(i + n_i, j + n_j)
                    P_neighbors = [P_neighbors; [i + n_i, j + n_j]];
                end
            end
        end

        % One neighbor means it's an end point.
        if size(P_neighbors, 1) == 1
            endness = true;

        % Two neighbors means it's an end point iff
        % they are 4-adjacent.
        elseif size(P_neighbors, 1) == 2
            first_neighbor = P_neighbors(1, 1: 2);
            second_neighbor = P_neighbors(2, 1: 2);

            % 4-adjacency between the two neighbors.
            if (first_neighbor(1) == second_neighbor(1) && ...
                abs(first_neighbor(2) - second_neighbor(2)) == 1) || ...
                (first_neighbor(2) == second_neighbor(2) && ...
                abs(first_neighbor(1) - second_neighbor(1)) == 1)
                endness = true;
            else
                endness = false;
            end
            
        % Any other case means it's not an end point.
        else
            endness = false;
        end
    
    % In 3D, an end point can also live in the rim of a surface.
    % Here, we follow Chris Pudney's approach, year 1998.
    elseif numel(point) == 3
        
        % Point coordinates.
        i = point{1};
        j = point{2};
        k = point{3};
        
        % Must be an object point.
        if ~V(i, j, k)
            endness = false;
            return;
        end
        
        % Point's 26-neighborhood.
        n_26 = V(i - 1: i + 1, j - 1: j + 1, k - 1: k + 1);
        
        % We define the 9 planes.
        planes(:, :, :, 1) = cat(3, [0 1 0; 0 1 0; 0 1 0], ...
                                    [0 1 0; 0 0 0; 0 1 0], ...
                                    [0 1 0; 0 1 0; 0 1 0]);
                                
        planes(:, :, :, 2) = cat(3, [0 0 0; 1 1 1; 0 0 0], ...
                                    [0 0 0; 1 0 1; 0 0 0], ...
                                    [0 0 0; 1 1 1; 0 0 0]);
                                
        planes(:, :, :, 3) = cat(3, [0 0 0; 0 0 0; 0 0 0], ...
                                    [1 1 1; 1 0 1; 1 1 1], ...
                                    [0 0 0; 0 0 0; 0 0 0]);
                                
        planes(:, :, :, 4) = cat(3, [1 0 0; 0 1 0; 0 0 1], ...
                                    [1 0 0; 0 0 0; 0 0 1], ...
                                    [1 0 0; 0 1 0; 0 0 1]);
                                
        planes(:, :, :, 5) = cat(3, [0 0 1; 0 1 0; 1 0 0], ...
                                    [0 0 1; 0 0 0; 1 0 0], ...
                                    [0 0 1; 0 1 0; 1 0 0]);
                                
        planes(:, :, :, 6) = cat(3, [1 1 1; 0 0 0; 0 0 0], ...
                                    [0 0 0; 1 0 1; 0 0 0], ...
                                    [0 0 0; 0 0 0; 1 1 1]);
                                
        planes(:, :, :, 7) = cat(3, [0 0 0; 0 0 0; 1 1 1], ...
                                    [0 0 0; 1 0 1; 0 0 0], ...
                                    [1 1 1; 0 0 0; 0 0 0]);
                                
        planes(:, :, :, 8) = cat(3, [1 0 0; 1 0 0; 1 0 0], ...
                                    [0 1 0; 0 0 0; 0 1 0], ...
                                    [0 0 1; 0 0 1; 0 0 1]);
                                
        planes(:, :, :, 9) = cat(3, [0 0 1; 0 0 1; 0 0 1], ...
                                    [0 1 0; 0 0 0; 0 1 0], ...
                                    [1 0 0; 1 0 0; 1 0 0]);
        
        % We check whether current point has exactly one neighbor in
        % any of these planes. If it does, it is an end point.
        for i = 1: numel(planes(1, 1, 1, :))
            if nnz(planes(:, :, :, i) .* n_26) == 1
                endness = true;
                return
            end
        end
        
        % If we get here, this point is not an end point.
        endness = false;
    
    % My brain won't handle 4D.
    else
        endness = 'naked lunch';
    end
end