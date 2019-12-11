% Our goal was programming a 3D thinning algorithm based upon
% MATLAB's morphological operations, replicating (kind of) what's
% done here: http://homepages.inf.ed.ac.uk/rbf/HIPR2/thin.htm
% input:
% V = 3D matrix (such as read_binvox.m's output).
% iters = iteration limit (optional)
% output: 3D matrix representing input's skeleton.
function skel = thinning3D(V, iters)

    % We ended up programming Palágyi's curve skel algorithm.
    % Reference: A parallel 3D 12-subiteration thinning algorithm.
    % Palágyi & Kuba,
    % Graphical Models and Image Processing, 1999, Elsevier.
    
    % 3D structuring elements for directional thinning.
    se = [];
    
    % Element count for each set.
    elem_count = [];
    
    % Format: (x, y, z, TX, element_number).
    % 1. T1.
    elem_count(1) = 16;
    se(:, :, 1, 1, 1) = [
        -1 -1 -1
         1  0  0
         0  0  0
    ];
    se(:, :, 2, 1, 1) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 1) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    % 2. T1 again.
    se(:, :, 1, 1, 2) = [
        -1 -1 -1
         0  1  0
         0  0  0
    ];
    se(:, :, 2, 1, 2) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 2) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    % 3. T1 again.
    se(:, :, 1, 1, 3) = [
        -1 -1 -1
         0  0  1
         0  0  0
    ];
    se(:, :, 2, 1, 3) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 3) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    % 4. T1 again.
    se(:, :, 1, 1, 4) = [
        -1 -1 -1
         0  0  0
         1  0  0
    ];
    se(:, :, 2, 1, 4) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 4) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    % 5. T1 again.
    se(:, :, 1, 1, 5) = [
        -1 -1 -1
         0  0  0
         0  1  0
    ];
    se(:, :, 2, 1, 5) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 5) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    % 6. T1 again.
    se(:, :, 1, 1, 6) = [
        -1 -1 -1
         0  0  0
         0  0  1
    ];
    se(:, :, 2, 1, 6) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 6) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    % 7. T1 again.
    se(:, :, 1, 1, 7) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    se(:, :, 2, 1, 7) = [
        -1 -1 -1
         1  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 7) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    % 8. T1 again.
    se(:, :, 1, 1, 8) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    se(:, :, 2, 1, 8) = [
        -1 -1 -1
         0  1  1
         0  1  0
    ];
    se(:, :, 3, 1, 8) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    % 9. T1 again.
    se(:, :, 1, 1, 9) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    se(:, :, 2, 1, 9) = [
        -1 -1 -1
         0  1  0
         1  1  0
    ];
    se(:, :, 3, 1, 9) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    % 10. T1 again.
    se(:, :, 1, 1, 10) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    se(:, :, 2, 1, 10) = [
        -1 -1 -1
         0  1  0
         0  1  1
    ];
    se(:, :, 3, 1, 10) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    % 11. T1 again.
    se(:, :, 1, 1, 11) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    se(:, :, 2, 1, 11) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 11) = [
        -1 -1 -1
         1  0  0
         0  0  0
    ];
    % 12. T1 again.
    se(:, :, 1, 1, 12) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    se(:, :, 2, 1, 12) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 12) = [
        -1 -1 -1
         0  1  0
         0  0  0
    ];
    % 13. T1 again.
    se(:, :, 1, 1, 13) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    se(:, :, 2, 1, 13) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 13) = [
        -1 -1 -1
         0  0  1
         0  0  0
    ];
    % 14. T1 again.
    se(:, :, 1, 1, 14) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    se(:, :, 2, 1, 14) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 14) = [
        -1 -1 -1
         0  0  0
         1  0  0
    ];
    % 15. T1 again.
    se(:, :, 1, 1, 15) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    se(:, :, 2, 1, 15) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 15) = [
        -1 -1 -1
         0  0  0
         0  1  0
    ];
    % 16. T1 again.
    se(:, :, 1, 1, 16) = [
        -1 -1 -1
         0  0  0
         0  0  0
    ];
    se(:, :, 2, 1, 16) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 1, 16) = [
        -1 -1 -1
         0  0  0
         0  0  1
    ];
    % 17. T2.
    elem_count(2) = 16;
    se(:, :, 1, 2, 1) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 1) = [
         1  0  0
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 1) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 18. T2 (again).
    se(:, :, 1, 2, 2) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 2) = [
         0  1  0
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 2) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 19. T2 (again).
    se(:, :, 1, 2, 3) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 3) = [
         0  0  1
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 3) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 20. T2 (again).
    se(:, :, 1, 2, 4) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 4) = [
         0  0  0
         1  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 4) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 21. T2 (again).
    se(:, :, 1, 2, 5) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 5) = [
         0  0  0
         0  1  1
         0  0  0
    ];
    se(:, :, 3, 2, 5) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 22. T2 (again).
    se(:, :, 1, 2, 6) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 6) = [
         0  0  0
         0  1  0
         1  0  0
    ];
    se(:, :, 3, 2, 6) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 23. T2 (again).
    se(:, :, 1, 2, 7) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 7) = [
         0  0  0
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 2, 7) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 24. T2 (again).
    se(:, :, 1, 2, 8) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 8) = [
         0  0  0
         0  1  0
         0  0  1
    ];
    se(:, :, 3, 2, 8) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 25. T2 (again).
    se(:, :, 1, 2, 9) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 9) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 9) = [
         1  0  0
         0  1  0
         0  0  0
    ];
    % 26. T2 (again).
    se(:, :, 1, 2, 10) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 10) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 10) = [
         0  1  0
         0  1  0
         0  0  0
    ];
    % 27. T2 (again).
    se(:, :, 1, 2, 11) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 11) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 11) = [
         0  0  1
         0  1  0
         0  0  0
    ];
    % 28. T2 (again).
    se(:, :, 1, 2, 12) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 12) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 12) = [
         0  0  0
         1  1  0
         0  0  0
    ];
    % 29. T2 (again).
    se(:, :, 1, 2, 13) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 13) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 13) = [
         0  0  0
         0  1  1
         0  0  0
    ];
    % 30. T2 (again).
    se(:, :, 1, 2, 14) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 14) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 14) = [
         0  0  0
         0  1  0
         1  0  0
    ];
    % 31. T2 (again).
    se(:, :, 1, 2, 15) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 15) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 15) = [
         0  0  0
         0  1  0
         0  1  0
    ];
    % 32. T2 (again).
    se(:, :, 1, 2, 16) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 2, 16) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 2, 16) = [
         0  0  0
         0  1  0
         0  0  1
    ];
    % 33. T3.
    elem_count(3) = 8;
    se(:, :, 1, 3, 1) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 3, 1) = [
        -1 -1 -1
         1  1  0
         0  0  0
    ];
    se(:, :, 3, 3, 1) = [
        -1 -1 -1
         0  0  0
         0  1  0
    ];
    % 34. T3 again.
    se(:, :, 1, 3, 2) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 3, 2) = [
        -1 -1 -1
         0  1  1
         0  0  0
    ];
    se(:, :, 3, 3, 2) = [
        -1 -1 -1
         0  0  0
         0  1  0
    ];
    % 35. T3 again.
    se(:, :, 1, 3, 3) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 3, 3) = [
        -1 -1 -1
         0  1  0
         1  0  0
    ];
    se(:, :, 3, 3, 3) = [
        -1 -1 -1
         0  0  0
         0  1  0
    ];
    % 36. T3 again.
    se(:, :, 1, 3, 4) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 3, 4) = [
        -1 -1 -1
         0  1  0
         0  0  1
    ];
    se(:, :, 3, 3, 4) = [
        -1 -1 -1
         0  0  0
         0  1  0
    ];
    % 37. T3 again.
    se(:, :, 1, 3, 5) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 3, 5) = [
        -1 -1 -1
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 3, 5) = [
        -1 -1 -1
         1  0  0
         0  1  0
    ];
    % 38. T3 again.
    se(:, :, 1, 3, 6) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 3, 6) = [
        -1 -1 -1
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 3, 6) = [
        -1 -1 -1
         0  0  1
         0  1  0
    ];
    % 39. T3 again.
    se(:, :, 1, 3, 7) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 3, 7) = [
        -1 -1 -1
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 3, 7) = [
        -1 -1 -1
         0  0  0
         1  1  0
    ];
    % 40. T3 again.
    se(:, :, 1, 3, 8) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 3, 8) = [
        -1 -1 -1
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 3, 8) = [
        -1 -1 -1
         0  0  0
         0  1  1
    ];
    % 41. T4.
    elem_count(4) = 3;
    se(:, :, 1, 4, 1) = [
        -1 -1 -1
        -1 -1 -1
         0  0  0
    ];
    se(:, :, 2, 4, 1) = [
         0 -1  0
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 4, 1) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 42. T4 again (v switch).
    se(:, :, 1, 4, 2) = [
        -1 -1 -1
         0 -1 -1
         0  0  0
    ];
    se(:, :, 2, 4, 2) = [
        -1 -1  0
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 4, 2) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 43. T4 again (v and w switch).
    se(:, :, 1, 4, 3) = [
        -1 -1 -1
         0 -1  0
         0  0  0
    ];
    se(:, :, 2, 4, 3) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 4, 3) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 44. T5.
    elem_count(5) = 4;
    se(:, :, 1, 5, 1) = [
        -1 -1  1
        -1 -1 -1
         0  0  0
    ];
    se(:, :, 2, 5, 1) = [
         0 -1  1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 5, 1) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 45. T5 again (v switch).
    se(:, :, 1, 5, 2) = [
        -1 -1  1
         0 -1 -1
         0  0  0
    ];
    se(:, :, 2, 5, 2) = [
        -1 -1  1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 5, 2) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 46. T5 again (z switch).
    se(:, :, 1, 5, 3) = [
        -1 -1  1
        -1 -1  1
         0  0  0
    ];
    se(:, :, 2, 5, 3) = [
         0 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 5, 3) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 47. T5 again (v and z switch).
    se(:, :, 1, 5, 4) = [
        -1 -1  1
         0 -1  1
         0  0  0
    ];
    se(:, :, 2, 5, 4) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 5, 4) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 48. T6.
    elem_count(6) = 4;
    se(:, :, 1, 6, 1) = [
         1 -1 -1
        -1 -1 -1
         0  0  0
    ];
    se(:, :, 2, 6, 1) = [
         1 -1  0
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 6, 1) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 49. T6 again (z switch).
    se(:, :, 1, 6, 2) = [
         1 -1 -1
         1 -1 -1
         0  0  0
    ];
    se(:, :, 2, 6, 2) = [
        -1 -1  0
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 6, 2) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 50. T6 again (v switch).
    se(:, :, 1, 6, 3) = [
         1 -1 -1
        -1 -1  0
         0  0  0
    ];
    se(:, :, 2, 6, 3) = [
         1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 6, 3) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 51. T6 again (z and v switch).
    se(:, :, 1, 6, 4) = [
         1 -1 -1
         1 -1  0
         0  0  0
    ];
    se(:, :, 2, 6, 4) = [
        -1 -1 -1
         0  1  0
         0  1  0
    ];
    se(:, :, 3, 6, 4) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 52. T7.
    elem_count(7) = 2;
    se(:, :, 1, 7, 1) = [
        -1 -1  0
        -1 -1  0
         0  0  0
    ];
    se(:, :, 2, 7, 1) = [
         0 -1  0
         0  1  1
         0  1  0
    ];
    se(:, :, 3, 7, 1) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 53. T7 (v switch).
    se(:, :, 1, 7, 2) = [
        -1 -1  0
         0 -1  0
         0  0  0
    ];
    se(:, :, 2, 7, 2) = [
        -1 -1  0
         0  1  1
         0  1  0
    ];
    se(:, :, 3, 7, 2) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 54. T8.
    elem_count(8) = 2;
    se(:, :, 1, 8, 1) = [
         0 -1 -1
         0 -1 -1
         0  0  0
    ];
    se(:, :, 2, 8, 1) = [
         0 -1  0
         1  1  0
         0  1  0
    ];
    se(:, :, 3, 8, 1) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 55. T8 (v switch).
    se(:, :, 1, 8, 2) = [
         0 -1 -1
         0 -1  0
         0  0  0
    ];
    se(:, :, 2, 8, 2) = [
         0 -1 -1
         1  1  0
         0  1  0
    ];
    se(:, :, 3, 8, 2) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 56. T9.
    elem_count(9) = 2;
    se(:, :, 1, 9, 1) = [
         1 -1  0
        -1 -1  0
         0  0  0
    ];
    se(:, :, 2, 9, 1) = [
         1 -1  0
         0  1  1
         0  1  0
    ];
    se(:, :, 3, 9, 1) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 57. T9 (z switch).
    se(:, :, 1, 9, 2) = [
         1 -1  0
         1 -1  0
         0  0  0
    ];
    se(:, :, 2, 9, 2) = [
        -1 -1  0
         0  1  1
         0  1  0
    ];
    se(:, :, 3, 9, 2) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 58. T10.
    elem_count(10) = 2;
    se(:, :, 1, 10, 1) = [
         0 -1  1
         0 -1 -1
         0  0  0
    ];
    se(:, :, 2, 10, 1) = [
         0 -1  1
         1  1  0
         0  1  0
    ];
    se(:, :, 3, 10, 1) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 59. T10 (z switch).
    se(:, :, 1, 10, 2) = [
         0 -1  1
         0 -1  1
         0  0  0
    ];
    se(:, :, 2, 10, 2) = [
         0 -1 -1
         1  1  0
         0  1  0
    ];
    se(:, :, 3, 10, 2) = [
         0  0  0
         0  1  0
         0  0  0
    ];
    % 60. T11.
    elem_count(11) = 1;
    se(:, :, 1, 11, 1) = [
        -1 -1 -1
        -1 -1  0
        -1 -1  0
    ];
    se(:, :, 2, 11, 1) = [
        -1 -1 -1
         0  1  0
         0  0  1
    ];
    se(:, :, 3, 11, 1) = [
        -1 -1 -1
         0  0  0
         0  1  0
    ];
    % 61. T12.
    elem_count(12) = 1;
    se(:, :, 1, 12, 1) = [
        -1 -1 -1
         0 -1 -1
         0 -1 -1
    ];
    se(:, :, 2, 12, 1) = [
        -1 -1 -1
         0  1  0
         1  0  0
    ];
    se(:, :, 3, 12, 1) = [
        -1 -1 -1
         0  0  0
         0  1  0
    ];
    % 62. T13.
    elem_count(13) = 1;
    se(:, :, 1, 13, 1) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 13, 1) = [
        -1 -1  0
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 13, 1) = [
        -1 -1  0
         0  0  1
         0  1  0
    ];
    % 63. T14.
    elem_count(14) = 1;
    se(:, :, 1, 14, 1) = [
        -1 -1 -1
        -1 -1 -1
        -1 -1 -1
    ];
    se(:, :, 2, 14, 1) = [
         0 -1 -1
         0  1  0
         0  0  0
    ];
    se(:, :, 3, 14, 1) = [
         0 -1 -1
         1  0  0
         0  1  0
    ];
    % Structuring element count.
    se_count = numel(se(1, 1, 1, :, 1));
    % Initially, aplying previous matrices as structuring elemnts
    % for hit-and-miss operation only remove U and S voxels.
    % Structuring elements must be rotated and hit-and-missed
    % successively through all 12 sub-iterations.
    % Thus, we must define a cyclic rotation order.
    % Check rot90mat.m.
    % US -> NE.
    rotations{1} = [3 2 2];
    % NE -> WD.
    rotations{2} = [1 1 1 2 2];
    % WD -> ES.
    rotations{3} = [1 2 2];
    % ES -> UW. 
    rotations{4} = [3 2 2 2];
    % UW -> ND.
    rotations{5} = [3 2 2 2];
    % ND -> SW.
    rotations{6} = [3 2 2];
    % SW -> UN.
    rotations{7} = [3 3 3 2 2];
    % UN -> ED.
    rotations{8} = [2 3 3];
    % ED -> NW.
    rotations{9} = [1 2];
    % NW -> UE.
    rotations{10} = [1 3 3 3]; 
    % UE -> SD.
    rotations{11} = [2 1 1];
    % SD -> US.
    rotations{12} = [1 1 1];
    % Rotations count (should be 12).
    rot_count = numel(rotations);
    skel = V;
    % Index of current structuring elements rotation.
    rot_index = 1;
    % Initial skeleton.
    prev_skel = skel;
    
    % Iteration limit is optional.
    % 1 iteration is useful for border detection.
    current_iteration = 0;
    w = 1;
    while true
        fprintf('rotation %d in iteration %d\n', w, ceil(w/12));
        w = w + 1;
        % Loop through all elements.
        compound_skel = zeros(size(skel));
        for i = 1: se_count
            % For most elements, we must hit and miss more than once.
            % To achieve this, we first compound the results and then
            % substract them from the partially calculated skeleton.
            compound_mask = zeros(size(skel));
            for j = 1: elem_count(i)
                compound_mask = compound_mask | bwhitmiss(skel, se(:, :, :, i, j));
            end
            compound_skel = compound_skel | compound_mask;
        end
        skel = skel - compound_skel;
        % Rotate all elements according to current rotation index.
        for i = 1: numel(rotations{rot_index})
            r = rotations{rot_index}(i);
            for j = 1: se_count
                for k = 1: elem_count(j)
                    se(:, :, :, j, k) = rot90mat(se(:, :, :, j, k), r);
                end
            end
        end
        % If we've reached the end of our rotations list, we must check
        % whether the skeleton went through unchanged. If so, we stop.
        % Otherwise, we start the process over.
        if rot_index == rot_count
            voxels_removed = nnz(prev_skel) - nnz(skel);
            fprintf('%d voxels removed\n', voxels_removed);
            if voxels_removed == 0
                break;
            else
                prev_skel = skel;
                rot_index = 1;
            end
            % Iterations limit.
            if nargin > 1
                current_iteration = current_iteration + 1;
                if current_iteration == iters
                    break
                end
            end            
        else
            rot_index = rot_index + 1;
        end
    end
end