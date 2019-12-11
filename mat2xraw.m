% Outputs an XRAW file from a MATLAB 3D binary matrix.
% XRAW specification:
% https://voxel.codeplex.com/wikipage?title=XRAW%20Format&version=7
% input: 3D binary matrix
% output: none, but outpurs matrix to given file_path.
function mat2xraw(V, file_path)

    % Create the .xraw file.
    fd = fopen(file_path, 'w');
    
    % Magic number.
    fwrite(fd, 'XRAW');
    
    % Channel data.
    fwrite(fd, 0);
    fwrite(fd, 4);
    fwrite(fd, 8);
    
    % Bits per index.
    bpi = 8;
    fwrite(fd, bpi);
    
    % Volume size, in little endian.
    vol_bytes = dec2hex(swapbytes(uint32(size(V, 3))), 8);
    fwrite(fd, hex2dec(vol_bytes(1:2)));
    fwrite(fd, hex2dec(vol_bytes(3:4)));
    fwrite(fd, hex2dec(vol_bytes(5:6)));
    fwrite(fd, hex2dec(vol_bytes(7:8)));
    vol_bytes = dec2hex(swapbytes(uint32(size(V, 2))), 8);
    fwrite(fd, hex2dec(vol_bytes(1:2)));
    fwrite(fd, hex2dec(vol_bytes(3:4)));
    fwrite(fd, hex2dec(vol_bytes(5:6)));
    fwrite(fd, hex2dec(vol_bytes(7:8)));
    vol_bytes = dec2hex(swapbytes(uint32(size(V, 1))), 8);
    fwrite(fd, hex2dec(vol_bytes(1:2)));
    fwrite(fd, hex2dec(vol_bytes(3:4)));
    fwrite(fd, hex2dec(vol_bytes(5:6)));
    fwrite(fd, hex2dec(vol_bytes(7:8)));
    
    % # Pallete colors (little endian 256).
    fwrite(fd, [0 1 0 0]);
    
    % Output voxel data.
    buffer_size = 1048576;
    buffer = zeros(1, buffer_size);
    buffer_index = 0;
    for i = 1: size(V, 1)
        for j = 1: size(V, 2)
            for k = 1: size(V, 3)
                buffer_index = buffer_index + 1;
                buffer(buffer_index) = V(i, j, k);
                if buffer_index == buffer_size
                    fwrite(fd, buffer);
                    buffer_index = 0;
                end
            end
        end
    end
    
    % Output last bytes.
    fwrite(fd, buffer(1: buffer_index));
    
    % Palette buffer data, only white.
    for i = 1: bpi * 256
        fwrite(fd, 255);
        fwrite(fd, 255);
        fwrite(fd, 255);
        fwrite(fd, 255);
    end
    
    % Close the file.
    fclose(fd);
end
