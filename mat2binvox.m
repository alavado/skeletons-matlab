% input: 3D matrix, name of output binvox file.
% output: Path of created binvox file.
function path = mat2binvox(V, fname)

    if(nargin < 2)
        fname = 'a.binvox';
    end
    
    % Binvox version.
    version = 1;
    
    % Check squareness of the matrix
    matrix_size = size(V);
    if size(matrix_size) ~= 3;
        error('Error: matrix must be 3D.');
    end
    
    % Grid dimensions.
    dim_1 = matrix_size(1);
    dim_2 = matrix_size(2);
    dim_3 = matrix_size(3);
    
    % Create the .binvox file.
    file_id = fopen(fname, 'w');

    % Header output.
    fprintf(file_id, '#binvox %d\n', version);
    
    % Grid dimensions.
    fprintf(file_id, 'dim %d %d %d\n', dim_1, dim_2, dim_3);
    
    % Header end keyword. We have no information about
    % original mesh, so no translation or scale are specified.
    fprintf(file_id, 'data\n');
    
    % Voxel compression.
    current_value = -1;
    count = 0;
    for i = 1: dim_1
        for j = 1: dim_2
            for k = 1: dim_3
                
                % Read current matrix index
                value = V(i, j, k);
                
                % If current_value differs from previous one,
                % we print to file and start counting again.
                % Also if we can't take it anymore.
                if current_value == -1
                    current_value = value;
                    count = 1;
                elseif value ~= current_value || count == 255
                    fwrite(file_id, current_value);
                    fwrite(file_id, count);
                    current_value = value;
                    count = 1;
                else
                    count = count + 1;
                end
            end
        end
    end
    
    % Last values.
    fwrite(file_id, current_value);
    fwrite(file_id, count);    
    fclose(file_id);
end