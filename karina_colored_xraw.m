function karina_colored_xraw(skel, output_file)
% Skeleton to XRAW with nice colors and stuff.

    % Keep only the largest connected component.
    CC = bwconncomp(skel);
    numOfPixels = cellfun(@numel,CC.PixelIdxList);
    [~, indexOfMax] = max(numOfPixels);
    skel = zeros(size(skel));
    skel(CC.PixelIdxList{indexOfMax}) = 1;

    % Paint the skeleton.
    skel = paint_skeleton(skel);
    
    skel(35, 406, 533) = 4;

    % Create the .xraw file.
    fd = fopen(output_file, 'w');
    
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
    for i = 3: -1: 1
        vol_bytes = dec2hex(swapbytes(uint32(size(skel, i))), 8);
        for j = 0: 3
            fwrite(fd, hex2dec(vol_bytes(2 * j + 1: 2 * j + 2)));
        end
    end
    
    % # Pallete colors (little endian 256).
    fwrite(fd, [0 1 0 0]);
    
    % Output voxel data.
    buffer_size = 1048576;
    buffer = zeros(1, buffer_size);
    buffer_index = 0;
    for i = 1: size(skel, 1)
        for j = 1: size(skel, 2)
            for k = 1: size(skel, 3)
                buffer_index = buffer_index + 1;
                buffer(buffer_index) = skel(i, j, k);
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
    for i = 1: 128
        fwrite(fd, [0 0 0 0]);
        
        % Normal voxel color.
        fwrite(fd, [255 255 255 255]);
        
        % End points color.
        fwrite(fd, [255 0 0 255]);
        
        % Joint points color.
        fwrite(fd, [255 0 0 255]);
        
        % ROI color.
        fwrite(fd, [255 255 0 255]);
    end
    
    % Close the file.
    fclose(fd);
end
