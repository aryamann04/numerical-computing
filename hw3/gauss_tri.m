function [l, u_diag, u_sup] = tridiag_LU(a, d, c)
%   a     : subdiagonal vector (length n-1) [A(i,i-1) for i=2:n]
%   d     : main diagonal vector (length n)    [A(i,i)   for i=1:n]
%   c     : superdiagonal vector (length n-1)   [A(i,i+1) for i=1:n-1]

    n = length(d);
    l = zeros(n-1, 1);   % multipliers for L (subdiagonal)
    u_diag = zeros(n, 1);% main diagonal of U
    u_sup = zeros(n-1, 1); % superdiagonal of U
    
    % initialization: first row of U is just the first row of A
    u_diag(1) = d(1);
    u_sup(1) = c(1);
    
    % for i = 2 to n, compute the multiplier and the new pivot.
    for i = 2:n
        % compute multiplier via L(i,i-1) = a(i-1) / u_diag(i-1)
        l(i-1) = a(i-1) / u_diag(i-1);
        
        % compute the main diagonal entry for U at row i via
        % u_diag(i) = d(i) - L(i,i-1)*c(i-1)
        u_diag(i) = d(i) - l(i-1) * c(i-1);
        
        % if not in the last row, set the superdiagonal entry from A
        if i < n
            u_sup(i) = c(i);
        end
    end
end

function y = forwardsub_tri(l, b)

    n = length(b);
    y = zeros(n,1);
    y(1) = b(1);
    
    % forward substitution for i=2 to n
    for i = 2:n
        y(i) = b(i) - l(i-1) * y(i-1);
    end
end

function x = backwardsub_tri(u_diag, u_sup, y)

    n = length(y);
    x = zeros(n,1);
    
    % back substitution, so start with the last row
    x(n) = y(n) / u_diag(n);
    
    for i = n-1:-1:1
        x(i) = (y(i) - u_sup(i) * x(i+1)) / u_diag(i);
    end
end

% generate a random tridiagonal matrix of order n
% a: subdiagonal (length n-1) for A(i,i-1)
% d: main diagonal (length n) for A(i,i)
% c: superdiagonal (length n-1) for A(i,i+1)
n = 10;  
a = randn(n-1,1);  
d = randn(n,1);
c = randn(n-1,1);

% d(1) = 1e-10;  % set A(1,1) to 1e-10

% construct the full tridiagonal matrix A explicitly
A = diag(d) + diag(c,1) + diag(a,-1);
b = randn(n,1);

% compute LU factorization for the tridiagonal matrix
[l, u_diag, u_sup] = tridiag_LU(a, d, c);

% solve L*y = b with forward substitution for unit lower bidiagonal L
y = forwardsub_tri(l, b);

% solve U*x = y with backward substitution for upper bidiagonal U
x = backwardsub_tri(u_diag, u_sup, y);
x1 = A\b; % MATLAB built-in for comparison

format long e
disp(x);
disp(x1);
norm(x - x1)