function x = truss_forces()
    alpha = sqrt(2)/2;
    A = zeros(13,13);
    b = zeros(13,1);

    % horizontal at A (1): alpha x1 - x4 - alpha x5 = 0
    A(1,1) =  alpha; 
    A(1,4) = -1;
    A(1,5) = -alpha;
    b(1)    =  0;

    % vertical at A (2): alpha x1 + x3 + alpha x5 = 0
    A(2,1) =  alpha;
    A(2,3) =  1;
    A(2,5) =  alpha;
    b(2)   =  0;

    % horizontal at B (3): x2 - x6 = 0
    A(3,2) =  1;
    A(3,6) = -1;
    b(3)   =  0;

    % vertical at B (4): x3 = 10
    A(4,3) =  1;
    b(4)   =  10;

    % horizontal at C (5): x4 - x8 = 0
    A(5,4) =  1;
    A(5,8) = -1;
    b(5)   =  0;

    % vertical at C (6): x7 = 0
    A(6,7) =  1;
    b(6)   =  0;

    % horizontal at D (7): x6 + alpha x5 - x10 - alpha x9 = 0
    A(7,6)  =  1;
    A(7,5)  =  alpha; 
    A(7,10) =  -1;
    A(7,9)  =  -alpha; 
    b(7)    =  0;

    % vertical at D (8): alpha x5 + alpha x9 = 15
    A(8,5) =  alpha;
    A(8,9) =  alpha;
    b(8)   =  15;

    % horizontal at E (9): x8 + alpha x9 - alpha x12 = 0
    A(9,8)  =  1;
    A(9,9)  =  alpha;
    A(9,12) =  -alpha;
    b(9)    =  0;

    % vertical at E (10): alpha x9 + x11 + alpha x12 = 0
    A(10,9)  =  alpha;
    A(10,11) =  1;
    A(10,12) =  alpha;
    b(10)    =  0;

    % horizontal at F (11): x10 - x13 = 0
    A(11,10) =  1;
    A(11,13) = -1;
    b(11)    =  0;

    % vertical at F (12): x11 = 10
    A(12,11) =  1;
    b(12)    =  10;

    % horizontal at G (13): x13 + alpha x12 = 0
    A(13,12) =  alpha;
    A(13,13) =  1;
    b(13)    =  0;

    % solve and display
    x = A \ b;
    disp(x);
end

function [A,x,b] = truss_forces_generalized(k, use_sparse)
    alpha = sqrt(2)/2; 
    
    n = 13 + (k - 1) * 8;  % number of bars and equations

    % set up sparse matrix for timing comparison
    if (use_sparse)
        A = sparse(n, n); 
        b = sparse(n, 1); 
    else
        A = zeros(n, n);
        b = zeros(n, 1);
    end

    % joint A_1 horizontal: alpha x1 - x4 - alpha x5 = 0
    A(1,1) =  alpha; 
    A(1,4) = -1;
    A(1,5) = -alpha;
    b(1)   =  0;

    % set up equations for each section
    for section = 1:k
        m = (section - 1) * 8;
        
        % joint A (vertical only)
        A(2 + m, 1 + m) = alpha;
        A(2 + m, 3 + m) = 1;
        A(2 + m, 5 + m) = alpha;
        b(2 + m) = 0;

        % joint B (horizontal and vertical)
        A(3 + m, 2 + m) = 1;
        A(3 + m, 6 + m) = -1;
        b(3 + m) = 0;

        A(4 + m, 3 + m) = 1;
        b(4 + m) = 10 + 20 * (section - 1); % load at B increases 10, 30, 50...

        % joint C (horizontal and vertical)
        A(5 + m, 4 + m) = 1;
        A(5 + m, 8 + m) = -1;
        b(5 + m) = 0;

        A(6 + m, 7 + m) = 1;
        b(6 + m) = 0;

        % joint D (horizontal and vertical)
        A(7 + m, 6 + m) = 1;
        A(7 + m, 5 + m) = alpha;
        A(7 + m, 10 + m) = -1;
        A(7 + m, 9 + m) = -alpha;
        b(7 + m) = 0;

        A(8 + m, 5 + m) = alpha;
        A(8 + m, 9 + m) = alpha;
        b(8 + m) = 20 + 20 * (section - 1); % load at D increases 20, 40, 60...

        % joint E (horizontal only, set for section = k separately)
        if section ~= k
            A(9 + m, 8 + m) = 1;
            A(9 + m, 9 + m) = alpha;
            A(9 + m, 12 + m) = -1;
            A(9 + m, 13 + m) = -alpha;
            b(9 + m) = 0;
        end
    end

    % final section E, F, G setup after loop
    m = (k - 1) * 8;

    % joint E (horizontal and vertical)
    A(9 + m, 8 + m) = 1;
    A(9 + m, 9 + m) = alpha;
    A(9 + m, 12 + m) = -alpha;
    b(9 + m) = 0;

    A(10 + m, 9 + m) = alpha;
    A(10 + m, 11 + m) = -1;
    A(10 + m, 12 + m) = alpha;
    b(10 + m) = 0;

    % joint F (horizontal and vertical)
    A(11 + m, 10 + m) = 1;
    A(11 + m, 13 + m) = -1;
    b(11 + m) = 0;

    A(12 + m, 11 + m) = 1;
    b(12 + m) = 10 + 20 * k;

    % end point G (horizontal force balance)
    A(13 + m, 12 + m) = alpha;
    A(13 + m, 13 + m) = 1;
    b(13 + m) = 0;

    % solve the system
    x = A \ b;
    disp(x);
end

function sparsity(A)

    [L, U, P] = lu(A);

    A_bandwidth = bandwidth(A); 
    A_inv_bandwidth = bandwidth(inv(A));
    L_bandwidth = bandwidth(L); 
    U_bandwidth = bandwidth(U); 

    fprintf('A bandwidth: %d\n', A_bandwidth); 
    fprintf('inv(A) bandwidth: %d\n', A_inv_bandwidth);
    fprintf('L bandwidth: %d\n', L_bandwidth);
    fprintf('U bandwidth: %d\n', U_bandwidth);

    % sparsity plot of A
    figure;
    subplot(2, 2, 1);
    spy(A);
    title('Sparsity Pattern of A');

    % sparsity plot of inverse of A
    A_inv = inv(A);
    subplot(2, 2, 2);
    spy(A_inv);
    title('Sparsity Pattern of A^{-1}');

    % sparsity plots of LU decomposition
    [L, U] = lu(A);
    subplot(2, 2, 3);
    spy(L);
    title('Sparsity Pattern of L');

    subplot(2, 2, 4);
    spy(U);
    title('Sparsity Pattern of U');

    % check if pivoting occurred
    if isequal(P, eye(length(A))) % if P is just the identity 
        disp('no pivoting occurred');
    else
        disp('pivoting occurred');
    end
end

function timing_comparison(A, b)

    % direct solution x = A \ b
    tic;
    x1 = A \ b;
    t1 = toc;
    fprintf('direct solve (A\\b): %.6f seconds\n', t1);

    % inverse method x = inv(A) * b
    tic;
    x2 = full(inv(A)) * b;
    t2 = toc;
    fprintf('inverse method (inv(A)*b): %.6f seconds\n', t2);

    % LU decomposition
    tic;
    [L, U] = lu(A);
    y = L \ b;
    x3 = U \ y;
    t3 = toc;
    fprintf('LU decomposition: %.6f seconds\n', t3);

    % sparse vs full comparison: plot running times for larger k values
    k_values = 1:30;
    times_full = zeros(size(k_values));
    times_sparse = zeros(size(k_values));

    for i = 1:length(k_values)

        % full matrix
        [A_full, ~, b_full] = truss_forces_generalized(i, false); 
        tic;
        A_full \ b_full;
        times_full(i) = toc;

        % sparse matrix
        [A_sparse, ~, b_sparse] = truss_forces_generalized(i, true);
        tic;
        A_sparse \ b_sparse;
        times_sparse(i) = toc;
    end
    
    % plot higher k values for sparse mode only
    k_values_sparse_only = 1:50; 
    times_sparse_only = zeros(size(k_values_sparse_only));
    
    for i = 1:length(k_values_sparse_only)
        [A_sparse, ~, b_sparse] = truss_forces_generalized(i, true);
        tic;
        A_sparse \ b_sparse;
        times_sparse_only(i) = toc;
    end

    figure;
    semilogy(k_values, times_full, '-o', 'DisplayName', 'Full Matrix');
    hold on;
    semilogy(k_values, times_sparse, '-o', 'DisplayName', 'Sparse Matrix');
    xlabel('k (Number of Sections)');
    ylabel('Execution Time (s)');
    title('Execution Time Comparison: Full vs Sparse');
    legend;
    grid on;

    figure;
    semilogy(k_values_sparse_only, times_sparse_only, '-o', ...
        'DisplayName', 'Sparse Matrix Only');
    xlabel('k (Number of Sections)');
    ylabel('Execution Time (s)');
    title('Execution Time for Sparse Matrix Only');
    legend;
    grid on;

end

% truss_forces()
[A,x,b] = truss_forces_generalized(10, false);
sparsity(A); 
% timing_comparison(A, b); 