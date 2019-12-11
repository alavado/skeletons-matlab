% Checks whether a point is simple in V.
% input:
%   V: N-dimensional binary matrix
%   point: cell array of integers (point coordinates)
% output: 1 if point is simple in V, 0 otherwise.
function simplicity = is_simple(V, point)

    % For 2D.
    % Euler characteristic of 3x3 neighboor graph with cycles of
    % length 3 removed must be equal to 1.
    if numel(point) == 2
        
        % Point coordinates.
        i = point{1};
        j = point{2};
        
        % Must be an object point.
        if ~V(i, j)
            simplicity = false;
            return;
        end
        
        % Each neighbor equals one vertex.
        vertices = zeros(1,8);
        vertices(1) = V(i - 1, j - 1);
        vertices(2) = V(i, j - 1);
        vertices(3) = V(i + 1, j - 1);
        vertices(4) = V(i + 1, j);
        vertices(5) = V(i + 1, j + 1);
        vertices(6) = V(i, j + 1);
        vertices(7) = V(i - 1, j + 1);
        vertices(8) = V(i - 1, j);
        
        % Inefficient but easily understandable adjacency matrix.
        edges = zeros(8, 8);
        
        % 4-adjacency edges.
        for n = 1: 7
            edges(n, n + 1) = vertices(n) * vertices(n + 1);   
        end
        edges(8, 1) = vertices(8) * vertices(1);
        
        % 8-adjacency edges and degenerate cycle removal.        
        if ~all(vertices(2: 4))
            edges(2, 4) = vertices(2) * vertices(4);
        end
        if ~all(vertices(4: 6))
            edges(4, 6) = vertices(4) * vertices(6);
        end
        if ~all(vertices(6: 8))
            edges(6, 8) = vertices(6) * vertices(8);
        end
        if ~all([vertices(1: 2) vertices(8)])
            edges(8, 2) = vertices(8) * vertices(2);
        end
        
        % Euler characteristic: |V| - |E|
        num_vertices = nnz(vertices);
        num_edges = nnz(edges);
        simplicity = (num_vertices - num_edges) == 1;
    
    % For 3D.
    % Bertrand & Malandain, 1993.
    % P is simple if the number of 26-connected components 26-adjacent
    % to it in its 26-neighborhood minus itself AND the number of
    % 6-connected components 6-adjacent to it in the background of
    % its 18-neighborhood minus itself are equal to 1.
    elseif numel(point) == 3
        
        % Must be an object point.
        if ~V(point{1}, point{2}, point{3})
            simplicity = false;
            return;
        end
        
        % All set.
        c = get_pudney_c(V, point);
        c_star = get_pudney_c_star(V, point);
        simplicity = c == 1 && c_star == 1;
    
    % Can't handle more than 3D.
    else
        simplicity = 'videodrome';
    end
end