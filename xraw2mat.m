% Parses an xraw voxel model to a MATLAB matrix, as
% outputted by MagicaVoxel 0.98.
% XRAW specification:
% https://voxel.codeplex.com/wikipage?title=XRAW%20Format&version=7
% input: xraw model path (.xraw file).
% output: N-D matrix, where N is the number of channels in the file.
function V = xraw2mat(file_path)

    % Open the .xraw file.
    fd = fopen(file_path, 'r');
    
    % Read some bytes.
    line = fgets(fd);
    buf_index = 1;
    buf = unicode2native(line);
    buf_size = numel(buf);
    
    % Processed bytes global index.
    byte_index = 0;
    
    % Header is always 24 bytes long.
    header_size = 24;
    header = zeros(1, header_size);
    header_read = false;
    while ~header_read
        while buf_index <= buf_size
            byte_index = byte_index + 1;
            if byte_index > header_size
                header_read = true;
                break;
            end
            header(byte_index) = buf(buf_index);
            buf_index = buf_index + 1;
        end
        if ~header_read
            line = fgets(fd);
            buf = unicode2native(line);
            buf_size = numel(buf);
            buf_index = 1;
        end
    end

    % "XRAW" ASCII characters.
    magic_number = header(1: 4);
    
    % Channels specs.
    % Channel data type. Always 0: unsigned integer?
    cdt = header(5);
    
    % Channel count. Always 4: ABGR?
    channels = header(6);
    
    % Bits per channel. Always 8?
    bpc = header(7);
    
    % Bits per index. Always 8?
    bpi = header(8);
    
    % Volume size.
    sizex = little_endian_hex_vec2big_endian_dec(header(9: 12));
    sizey = little_endian_hex_vec2big_endian_dec(header(13: 16));
    sizez = little_endian_hex_vec2big_endian_dec(header(17: 20));
    
    % Always 256 (0100)?
    palette_colors = little_endian_hex_vec2big_endian_dec(header(21: 24));
    
    % Voxel data.
    voxel_count = sizex * sizey * sizez;
    if bpi == 0
        voxel_data_size = voxel_count * bpc * channels / 8;
    else
        voxel_data_size = voxel_count * bpi / 8;
    end
    
    % We will reshape this array later on.
    V = zeros(1, voxel_data_size);
    vi = 1;
    voxel_data_read = false;
    while ~voxel_data_read
        while buf_index <= buf_size
            if byte_index > header_size + voxel_data_size
                voxel_data_read = true;
                break;
            else
                V(vi) = buf(buf_index);
                vi = vi + 1;
                buf_index = buf_index + 1;
            end
            byte_index = byte_index + 1;
        end
        if ~voxel_data_read
            line = fgets(fd);
            buf = unicode2native(line);
            buf_size = numel(buf);
            buf_index = 1;
        else
            break;
        end
    end
    
    % Final Matlab matrix conversion.
    V = reshape(V, sizex, sizey, sizez);
    
    return;
    
    % Palette data.
    palette_data_size = bpi * palette_colors;
    total_data_size = header_size + voxel_data_size + palette_data_size;
    palette_data_read = false;
    bytespal = 0;
    while ~palette_data_read
        for i = byte_index: bytes_read
            bytespal = bytespal + 1;
            byte_index = byte_index + 1;
            if byte_index > total_data_size
                palette_data_read = true;
                break;
            end
        end
        if ~palette_data_read
            line = fgets(fd);
            if line == -1
                break;
                error('XRAW file does not meet specs');
            end
            bytes = unicode2native(line);
            bytes_read = bytes_read + numel(bytes);
        end
    end
    
    fclose(fd);
end

