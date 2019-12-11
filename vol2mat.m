% This function parses a .VOL voxel model to a MATLAB matrix.
% VOL specification: http://www.tc18.org/code_data_set/Code/Old_releases/libvol-README
% input:
%   model_path; VOL model path (.VOL file, for a model of dimension X*Y*Z).
%   line_break: optional specification of line breaks in file.
% output:
%   3D matrix, size = [X Y Z].
function V = vol2mat(model_path, line_break)
    
    % This is the line break that came in Gabriella's files.
    if(nargin < 2)
        line_break = '.\r\n';
    end

    % Open the .vol  file.
    file_id = fopen(model_path, 'r');

    % VOL header is a list of lines containing metadata.
    % X dimension
    line = fgets(file_id);
    [coord, val] = strtok(line);
    x_dim = str2num(val);
    
    % Y dimension.
    line = fgets(file_id);
    [coord, val] = strtok(line);
    y_dim = str2num(val);
    
    % Z dimension.
    line = fgets(file_id);
    [coord, val] = strtok(line);
    z_dim = str2num(val);
    
    % The rest of the header is not relevant for this thesis.
    while ~strcmp(line, sprintf(line_break))
        line = fgets(file_id);
    end
    
    % Also ignore the mysterious line that came after the header.
    line = fgets(file_id);
    
    % Array for each byte.
    V_array = zeros(1, x_dim * y_dim * z_dim);
    V_array_index = 1;
    
    % Bytes currently read.
    read_bytes = 0;
    
    % Now comes the data bytes.
    while 1
        
        % Read some data.
        line = fgets(file_id);
        
        % Break if EOF.
        if line == -1
            break;
        end
        
        % Line to bytes.
        bytes = unicode2native(line);
        num_bytes = numel(bytes);
       	read_bytes = read_bytes + num_bytes;
        
        for i = 1 : num_bytes
            V_array(V_array_index) = bytes(i);
            V_array_index = V_array_index + 1;
        end
    end
    
    % Matrix dimension is off.
    if read_bytes ~= x_dim * y_dim * z_dim
        error('voxel count does not match specified header dimensions');
    end
    
    % Stack the array into a cubic model.
    V = [];
    area = y_dim * z_dim;
    for i = 1: x_dim
        slice = vec2mat(V_array((i - 1) * area + 1: i * area), x_dim);
        V = cat(3, slice, V);
    end
    
    % Fix for viewvox's default orientation.
    V = rot90mat(V, 2);
    for i = 1: x_dim
        V(:, :, i) = fliplr(V(:, :, i));
    end
    
    return;
    
    % This flag indicates whether current byte is a value byte
    % or a count byte, following binvox's compression method.
    % This is necessary as MATLAB cannot handle byte streams, so
    % we can be surprised with a new line at any moment.
    value_byte = true;
    while line ~= -1
        line = fgets(file_id);
        if line == -1
            break
        end
        
        % String to byte array casting.
        bytes = unicode2native(line);
        num_bytes = numel(bytes);
        
        % We read each byte in current line.
        for i = 1: num_bytes
            if value_byte
                
                % Values must be binary.
                if bytes(i) ~= 1 && bytes(i) ~= 0
                    error(['Error: invalid value <' bytes(i) '>.']);
                end
                current_value = bytes(i);
            else
                for j = 1: bytes(i)
                    grid(grid_index) = current_value;
                    grid_index = grid_index + 1;
                end
                read_voxels = read_voxels + uint32(bytes(i));
            end
            
            % Switch byte type flag.
            value_byte = ~value_byte;
        end
    end

    % We must have read each and every voxel.
    if read_voxels ~= grid_size
        error(['Error: missing voxels, found ' num2str(read_voxels)...
            ' expecting ' num2str(grid_size) '.']);
    end

    % Finally, we divide and stack the grid array to form the volume.
    V = [];
    area = depth * height;
    for i = 1: depth
        slice = vec2mat(grid((i - 1) * area + 1: i * area), depth);
        V = cat(3, slice, V);
    end
    
    % Fix this later maybe.
    V = rot90mat(V, 2);
    for i = 1: depth
        V(:, :, i) = fliplr(V(:, :, i));
    end

    % Bye bye.
    disp('binvox2mat: Ok.');
    fclose(file_id);
end