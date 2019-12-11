% Rotates a square matrix in 90 degrees.
% input: Some square matrix and rotation axis dimension (1, 2 or 3).
% output: Rotated matrix.
function V = rot90mat(M, axis)
    
    M_size = size(M);
    
    switch axis
        case 1
            V = permute(M, [3 2 1]);
            for i = 1: M_size(3)
                V(:, :, i) = flipud(V(:, :, i));
            end
        case 2
            V = permute(M, [1 3 2]);
            for i = 1: M_size(2)
                V(:, :, i) = fliplr(V(:, :, i));
            end
        case 3
            % This is already implemented.
            V = imrotate(M, 90);
        otherwise
            error('Error: invalid rotation axis.');
    end
end