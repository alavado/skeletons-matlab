% Multiplies z-scale by 3.
function rescaled_neuron = rescale_neuron(neuron, dimension, scaling_factor)
    
    % Read shortest dimension if no dimension is specified.
    [x, y, z] = size(neuron);
    if nargin < 2
        if x <= y && x <= z
            dimension = 1;
        elseif y <= x && y <= z
            dimension = 2;
        else
            dimension = 3;
        end
    end
    
    % Scale up in dimension.
    if dimension == 1
        rescaled_neuron = zeros(scaling_factor * x, y, z);
        for i = 0: x - 1
            rescaled_neuron(3 * i + 1, :, :) = neuron(i + 1, :, :);
            rescaled_neuron(3 * i + 2, :, :) = neuron(i + 1, :, :);
            rescaled_neuron(3 * i + 3, :, :) = neuron(i + 1, :, :);
        end
    elseif dimension == 2
        rescaled_neuron = zeros(x, scaling_factor * y, z);
        for i = 0: y - 1
            rescaled_neuron(:, 3 * i + 1, :) = neuron(:, i + 1, :);
            rescaled_neuron(:, 3 * i + 2, :) = neuron(:, i + 1, :);
            rescaled_neuron(:, 3 * i + 3, :) = neuron(:, i + 1, :);
        end
    else
        rescaled_neuron = zeros(x, y, scaling_factor * z);
        for i = 0: z - 1
            for j = 1: scaling_factor
                rescaled_neuron(:, :, scaling_factor * i + j) = neuron(:, :, i + 1);
            end
        end
    end
end

