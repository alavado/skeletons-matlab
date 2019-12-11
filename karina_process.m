% File paths list.
stacks = {'D:\RSI\Karina\Paper_PpO_light\Karina\tallarinolebia_stack_0-255_1.tif'};

% Open the results file.
f = fopen('exp.txt','w');
fprintf(f, 'Archivo,Largo total,Bifurcaciones,Puntas\n');

xyzSizeFactor = 0.104;
% Skeletonize all stacks.
for i = 1: numel(stacks)
    [m1, m2, m3, m4, m5, m6, m7, m8] = get_skel_metrics(stacks{i}, [0 0 0]);
    fprintf(f, '%s,%d,%d,%d\n', stacks{i}, m3 * xyzSizeFactor, m7, m8);
end