 function [A,q] = house(A)
 %
 % function [A,q] = house(A)
 %
 % Perform QR decomposition using Householder reflections
 % R is returned in the upper triangular of A
 % Q is not returned except implicitly via the reflection vectors:
 % Transformations are of the form Q_k = I - 2u_k(u_k^T) 
 % First entry of vector u_k is returned in q(k), remaining entries in 
 % A(k+1:m,k). Assume m>n.

    [m,n]=size(A); q = zeros(1,n);
    for k = 1:n
        % define u of length = m-k+1
        z = A(k:m,k);
        e1 = [1; zeros(m-k,1)];
        u = z+sign(z(1))*norm(z)*e1; u = u/norm(u);
        % update nonzero part of A by I-2uu^T
        A(k:m,k:n) = A(k:m,k:n)-2*u*(u'*A(k:m,k:n));
        % store u
        q(k) = u(1);
        A(k+1:m,k) = u(2:m-k+1);
    end
end

function [A,q] = house_mod_q4(A)

    [m,n]=size(A); q = zeros(1,n);

    for k = 1:n
        % define u of length = m-k+1
        z = A(k:m,k);
        e1 = [1; zeros(m-k,1)];
        u = z+sign(z(1))*norm(z)*e1; u = u/norm(u);

        % only multiply householder reflector k onto 
        % subcolumns k+1 through n 
        A(k:m,k+1:n) = A(k:m,k+1:n)-2*u*(u'*A(k:m,k+1:n));

        % explicitly store the diagonal value 
        A(k,k) = -sign(z(1))*norm(z); 

        % store u
        q(k) = u(1);
        A(k+1:m,k) = u(2:m-k+1);
    end
end

function [A,q] = house_mod_q5(A)

    [m,n]=size(A); q = zeros(1,n);
    for k = 1:n
        % define u of length = m-k+1
        z = A(k:m,k);
        e1 = [1; zeros(m-k,1)];
        u = z+sign(z(1))*norm(z)*e1; u = u/norm(u);
        % update nonzero part of A by I-2uu^T
        A(k:m,k:n) = A(k:m,k:n)-2*(u*u')*A(k:m,k:n);
        % store u
        q(k) = u(1);
        A(k+1:m,k) = u(2:m-k+1);
    end
end

function timing_comparison()
    nvec = 10:10:200; % from n=10 to n=200 in steps of 10
    original_times = zeros(size(nvec)); 
    modified_times = zeros(size(nvec));
    
    for n = 1:length(nvec)
        A = randn(nvec(n) * 10, nvec(n)); 
        % time original function
        tic; 
        [A1, q1] = house(A); 
        original_times(n) = toc; 
        % time modified function 
        tic; 
        [A2, q2] = house_mod_q5(A); 
        modified_times(n) = toc;
    end
    
    % standard plot
    figure;
    plot(nvec, original_times, 'r-', nvec, modified_times, 'b-');
    xlabel('n'); 
    ylabel('time (seconds)'); 
    title(['timing comparison between original and modified' ...
        ' householder functions'])
    legend('original', 'modified')

    % plot log-log 
    figure;
    loglog(nvec, original_times, 'r-', nvec, modified_times, 'b-');
    xlabel('n'); 
    ylabel('time (seconds)'); 
    title(['timing comparison between original and modified' ...
        ' householder functions (semi-log scale)'])
    legend('original', 'modified')

end

function [Q,R] = myqr(A, flag)
    
    [m,n] = size(A); 

    % get the householder reduction to isolate Rhat 
    [A_house, q] = house(A); 
    Rhat = triu(A_house(1:n,1:n));
    
    if nargin() == 2
        % reduced factorization
        Q = eye(m,n); 
        R = Rhat; 
    else
        % regular factorization 
        Q = eye(m); 
        R = [Rhat; zeros(m-n,n)]; % append zeroes
    end 

    % Q = Q(1) * ... * Q(n-1) * Q(n) * I, so start at n 
    for k=n:-1:1 
        % construct u_k by finding its first element (in q(k))
        % and the rest in A_house(k+1:m, k)
        u = [q(k); A_house(k+1:m,k)];

        % update Q efficiently by only updating rows k through m
        Q(k:m,:) = Q(k:m,:) - 2 * u * (u' * Q(k:m,:)); 
    end
end

function QR_timing_comparison()
    nvec = 10:10:200; % from n=10 to n=200 in steps of 10
    matlab_qr_times = zeros(size(nvec)); 
    myqr_times = zeros(size(nvec));
    
    % also initialize a vector to hold the ratio of times
    ratio = zeros(size(nvec)); 
    
    for n = 1:length(nvec)
        A = randn(nvec(n) * 10, nvec(n)); 
        % time Matlab's qr function
        tic; 
        [Q, R] = qr(A); 
        matlab_qr_times(n) = toc;

        % time myqr function 
        tic; 
        [my_Q, my_R] = myqr(A); 
        myqr_times(n) = toc;
        
        % store ratio 
        ratio(n) = myqr_times(n)/matlab_qr_times(n);
    end
    
    figure;
    yyaxis left
    plot(nvec, matlab_qr_times, 'r-', nvec, myqr_times, 'b-');
    ylabel('time (seconds)');

    yyaxis right
    plot(nvec, ratio, 'k-*');
    ylabel('ratio (myQR/Matlab QR)');

    xlabel('n');
    title('timing comparison: myQR vs Matlab QR');
    legend('Matlab QR', 'myQR', 'ratio');

    figure;
    yyaxis left
    loglog(nvec, matlab_qr_times, 'r-', nvec, myqr_times, 'b-');
    ylabel('time (seconds)');

    yyaxis right
    loglog(nvec, ratio, 'k-*');
    ylabel('ratio (myQR/Matlab QR)');

    xlabel('n');
    title('log-Log timing comparison: myQR vs Matlab QR');
    legend('Matlab QR', 'myQR', 'ratio');

    fprintf("avg ratio: %e", mean(ratio)); 
end

% testing myqr 

[Q,R] = myqr(A); 
resid_norm = norm(A - Q*R); 
ortho_norm = norm(Q'*Q - eye(size(Q,2))); 

[Qhat, Rhat] = myqr(A,0); 
resid_hat_norm = norm(A - Qhat*Rhat); 
ortho_hat_norm = norm(Qhat'*Qhat - eye(size(Qhat,2))); 

fprintf('regular QR: residual norm = %e, orthogonality norm = %e\n', ...
    resid_norm, ortho_norm);
fprintf('reduced QR: residual norm = %e, orthogonality norm = %e\n\n', ...
    resid_hat_norm, ortho_hat_norm);

fprintf('R is upper triangular: %s\n', mat2str(istriu(R)));
fprintf('Rhat is upper triangular: %s\n', mat2str(istriu(Rhat)));

QR_timing_comparison(); 

[Aout1, ~] = house(A); 
% [Aout2, ~] = house_mod_q4(A);
[Aout3, ~] = house_mod_q5(A);
norm(Aout1 - Aout3)
% timing_comparison()




