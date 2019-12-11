% Little endian hex to decimal conversion.
% There must be a better way to achieve this.
% input: bytes vector (example: [45 0 0 0]).
% output: decimal number (45).
function d = little_endian_hex_vec2big_endian_dec(v)
    v_as_str = '';
    for i = 1: numel(v)
        v_as_str = strcat(v_as_str, num2str(v(i), '%02x'));
    end
    d = swapbytes(uint32(hex2dec(v_as_str)));
end
