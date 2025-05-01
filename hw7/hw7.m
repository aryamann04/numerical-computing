% Aryaman Nagpal 
% Numerical Computing HW 6

%------------%
% Question 1 % 
%------------%

function [v, mu, resids] = power_method(A, v0)
    
    v = v0 / norm(v0); % normalize initial guess 
    resids = zeros(100, 1); % vector for residuals

    for k = 1:100 

        v = A * v; 
        v = v / norm(v); 
        mu = v' * A * v; 
        resids(k) = norm(A * v - mu * v); 
        
        fprintf("iteration %d:\n", k); 
        fprintf("mu_%d = %20.16e\n", k, mu)
        fprintf("eigenvalue residual norm: %20.16e\n\n", resids(k));  
        
        if (resids(k) < 1e-14)
            break; 
        end 

    end 
end

function [v, mu, resids] = inverse_iter(A, v0, alpha)
    
    v = v0 / norm(v0); % normalize initial guess 

    % compute LU decomposition outside loop for efficiency
    [L, U] = lu(A - alpha * eye(size(A))); 
    
    resids = zeros(100, 1); % vector for residuals

    for k = 1:100 
        
        % use LU factors to solve (A - alpha * I)v = v_k-1
        v = U \ (L \ v); 
        v = v / norm(v); 
        mu = v' * A * v; 
        resids(k) = norm(A * v - mu * v); 
        
        fprintf("iteration %d:\n", k); 
        fprintf("mu_%d = %20.16e\n", k, mu)
        fprintf("eigenvalue residual norm: %20.16e\n\n", resids(k)); 
        
        if (resids(k) < 1e-14)
            break; 
        end 
        
    end 
end

function [v, mu, resids] = rayleigh_quot_iter(A, v0)
    
    v = v0 / norm(v0); % normalize initial guess  
    mu = v' * A  * v; % initial rayleigh quotient 
    resids = zeros(100, 1); % vector for residuals

    for k = 1:100 
        
        % solve (A - mu_k-1 * I)v = v_k-1
        M = A - mu * eye(size(A)); 
        v_tilde = M \ v; 
        v = v_tilde / norm(v_tilde); 

        % update quotient 
        mu = v' * A * v; 
        resids(k) = norm(A * v - mu * v); 

        fprintf("iteration %d:\n", k); 
        fprintf("mu_%d = %20.16e\n", k, mu)
        fprintf("eigenvalue residual norm: %20.16e\n\n", resids(k)); 
        
        if (resids(k) < 1e-14)
            break; 
        end 
        
    end 
end

function T = tridiagExample(N,nonsym)
%
% set up the tridiagonal matrix for hw6
% when nonsym = 0, this is the matrix shown on AG p.94-95, but without the 1/h^2 factors
% when nonsym is not 0, this is a matrix related to that discussed in 
% AG, p.304-305, but without the 1/(2h) factors
%
T = sparse(N);
T(1,1) = 2;
T(1,2) = -1 - nonsym;
for k=2:N-1
    T(k,k-1) = -1 + nonsym;
    T(k,k) = 2;
    T(k,k+1) = -1 - nonsym;
end
T(N,N-1) = -1 + nonsym;
T(N,N) = 2;
end 

function q1_main()
    
    T = tridiagExample(100, 0.1); 
    v0 = ones(100, 1); 
    
    [~, ~, resids_power] = power_method(T, v0); 
    [~, ~, resids_inv] = inverse_iter(T, v0, 4);
    [~, ~, resids_rayleigh] = rayleigh_quot_iter(T, v0);

    figure;

    semilogy(1:100, resids_power, 'r', 'LineWidth', 2); hold on;
    semilogy(1:100, resids_inv, 'g', 'LineWidth', 2);
    semilogy(1:100, resids_rayleigh, 'b', 'LineWidth', 2);

    legend('power method', 'inverse iteration', 'rayleigh quotient iteration');
    xlabel('iteration');
    ylabel('residual norm (log scale)');
    title('residual norms vs iteration count');
    grid on;
end 


%------------%
% Question 2 % 
%------------%


function q2_main()
        
    n = 11; 

    A = rand(n); 
    evalues = eig(A); 
    lambda = real(evalues(1)); 

    b = rand(n, 1); 
    A_lambda = A - lambda * eye(n); 
    x = (A_lambda) \ b; 

    x = x / norm(x); 

    fprintf("eigenvalue-vector residual: %20.16e\n", norm(A * x - lambda * x)); 

    [~, S, V] = svd(A_lambda); 

    % find index i of minimum singular value, which corresponds to 
    % the eigenvalue lambda
    [~, i] = min(diag(S)); 

    % extract that column (from the right singular vectors V) at index i
    sv = V(:, i); 

    fprintf("minimum singular value: %20.16e\n", S(i, i)); 
    fprintf("svd eigenvector residual: %20.16e\n\n", norm(A * sv - lambda * sv));

end


%------------%
% Question 3 % 
%------------%


%------%
% 3(a) % 
%------%

function [Q,R] = qrPosDiagR(A)

    [Q, R] = qr(A, 0);
    [~, n] = size(A); 

    for i = 1:n
        if R(i,i) < 0
            R(i, :) = R(i, :) * -1; 
            Q(:, i) = Q(:, i) * -1; 
        end
    end
end

%------------%
% 3(b), 3(c) % 
%------------%

function [U_k, D_U, V_k] = svd_block_power_method(A, r, tol, print_details)

    [~, n] = size(A); 
    k = 0; 
    V_k = eye(n, r); 

    while (true)

        % incremement k
        k = k + 1; 

        U_hat = A * V_k; 
        [Q_U, R_U] = qrPosDiagR(U_hat); 
        U_k = Q_U; 

        V_hat = A' * U_k; 
        [Q_V, R_V] = qrPosDiagR(V_hat); 
        V_k = Q_V; 

        resid_U = norm(triu(R_U, 1)); 
        resid_V = norm(triu(R_V, 1)); 

        if (max(resid_U, resid_V) < tol)
            break; 
        end 
    end 

    % recover only the diagonal part of R_U to eliminate small
    % errors in off-diagonal entries (returned as Sigma)
    D_U = diag(diag(R_U));  

    if print_details ~= 0

        format longE;

        % display the diagonals of R_U, R_V
        fprintf("diag(R_U) = \n"); 
        disp(diag(R_U)); 
        fprintf("diag(R_V) = \n"); 
        disp(diag(R_V)); 
    
        % display the first (largest) r singular values of A
        [~, S, ~] = svd(A);
        fprintf("diag(S) = \n"); 
        disp(diag(S(1:r,1:r))); 
 
        % part (c): print out the diagonal residuals and check if both <= tol 

        % recover only the diagonal part of R_V to eliminate small
        % errors in off-diagonal entries 
        D_V = diag(diag(R_V)); 

        % ||A * V_k - U_k * D_U||
        svd_block_resid_U = norm(A * V_k - U_k * D_U); 
        fprintf("residual for U: %20.16e\n", svd_block_resid_U);

        % ||A' * U_k - V * D_V||
        svd_block_resid_V = norm(A' * U_k - V_k * D_V); 
        fprintf("residual for V: %20.16e\n", svd_block_resid_V);

        % check <= tol 
        if (svd_block_resid_U <= tol && svd_block_resid_V <= tol) 
            fprintf("both residuals less than or equal to tolerance\n\n"); 
        else 
            fprintf("failure: at least one residual greater than tolerance\n\n");
        end 
    end 

end


%------------%
% Question 4 % 
%------------%


%-------------%
% main for Q4 % 
%-------------%

function q4_main()

    A = imread('/Users/aryaman/Downloads/WallOfWindows.jpg', 'jpg'); 
    % image(A);

    A_red = double(A(:,:,1)); 
    A_green = double(A(:,:,2)); 
    A_blue = double(A(:,:,3)); 
    
    % part (a), generating images using MATLAB's svd

    % generate_images(A_red, A_green, A_blue, 0, 0);
    % fprintf("MATLAB's svd, completed"); 
     
    % part (b), generating images using block power method 
    
    generate_images(A_red, A_green, A_blue, 1, 1e5);
    fprintf("block power method, tol = 100000, completed\n"); 

    generate_images(A_red, A_green, A_blue, 1, 1e4);
    fprintf("block power method, tol = 10000, completed\n"); 

    generate_images(A_red, A_green, A_blue, 1, 1e3);
    fprintf("block power method, tol = 1000, completed\n"); 

    generate_images(A_red, A_green, A_blue, 1, 1e2);
    fprintf("block power method, tol = 100, completed\n");  
    
end

% helper function:  generates the images for ranks r = 10, 20, 30, 100. The
%                   argument 'alg' chooses whether the block power method 
%                   (1) or MATLAB svd method (0) is used for the 
%                   decomposition. The argument 'tol is a tolerance 
%                   parameter passed to the block power method (ignored if
%                   MATLAB's svd is used)

function generate_images(A_red, A_green, A_blue, alg, tol)
    
    figure; 
    ranks = [10, 20, 30, 100]; 
    
    for i = 1:4
        
        r = ranks(i); 

        % get the nearest r rank matrix using a helper function that 
        % internally computes the SVD decomposition of each matrix 
        
        fprintf("starting alg = %d, tol = %20.16e, r = %d\n", alg, tol, r);

        tic; 
        A_red_r = nearest_r_rank_matrix(A_red, r, alg, tol); 
        A_green_r = nearest_r_rank_matrix(A_green, r, alg, tol); 
        A_blue_r = nearest_r_rank_matrix(A_blue, r, alg, tol);
        time = toc; 

        fprintf("ended alg = %d, tol = %20.16e, r = %d, time = %20.16e\n\n", alg, tol, r, time); 

        % recombine 
        A_r = uint8(cat(3, A_red_r, A_green_r, A_blue_r)); 

        % display the image 
        subplot(2, 2, i);  
        image(A_r);
        s = sprintf('r = %d', ranks(i)); 
        
        % add tolerance level to title if using the block power method
        if alg == 1 
            s = strcat(s, sprintf(', tol = %f', tol));
        end 

        title(s); 
    end 

end

% helper function:  for part (a), takes in matrix A and returns the 
%                   nearest r rank matrix to A, internally computes 
%                   SVD decomposition. The method used to compute the
%                   SVD is decided by 'alg': if 0, MATLAB's svd is used,
%                   if 1, then the block power method is used instead 

function A = nearest_r_rank_matrix(A, r, alg, tol)
    
    % get the SVD decomposition of A:
    % if alg is 0, use MATLAB's svd, otherwise use block power method

    if alg == 0

        [U, S, V] = svd(A, 0); 

    elseif alg == 1
        
        % note: the singular value matrix S is chosen to be D_U, the 
        % diagonal of R_U
        [U, S, V] = svd_block_power_method(A, r, tol, 0);

    else

        error("nearest_r_rank_matrix: 'alg' must be 0 or 1");

    end 

    % construct the nearest r rank matrix
    A = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)'; 

end 

q4_main(); 


