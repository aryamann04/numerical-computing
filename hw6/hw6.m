% Aryaman Nagpal 
% Numerical Computing HW 6

% Q4

function svd_plot()
    % parametrize unit circle via polar coordinates:
    %   x = rcos(t), y = rsin(t) 
    %   r = 1 (unit circle), theta in [0, 2pi)
    
    t = linspace(0, 2 * pi); 
    x = cos(t); 
    y = sin(t); 
    
    % generate random 2x2 matrix 
    A = randn(2); 
    
    % multiply onto (x,y) and isolate x and y components
    A_xy = A * [x; y];  
    A_x = A_xy(1, :); 
    A_y = A_xy(2, :); 
    
    % perform SVD and isolate u1, u2, v1, v2, sigma_1, sigma_2
    [U,S,V] = svd(A); 
    u1 = U(:, 1); 
    u2 = U(:, 2); 
    v1 = V(:, 1); 
    v2 = V(:, 2); 
    sigma_1 = S(1,1); 
    sigma_2 = S(2,2); 
    
    % generate side-by-side plots of the unit circle and ellipse 
    
    % first construct the unit circle plot
    subplot(1,2,1); 
    plot(x, y);
    axis equal; 
    xlabel('x'); 
    xlabel('y'); 
    title('unit circle plot');
    
    % add v1, v2 vectors to unit circle plot 
    hold on; 
    plot([0, v1(1)], [0, v1(2)], 'LineWidth', 2);
    text(v1(1), v1(2), '  v1', 'FontSize', 12);
    plot([0, v2(1)], [0, v2(2)], 'LineWidth', 2);
    text(v2(1), v2(2), '  v2', 'FontSize', 12);
    hold off; 
    
    % now construct the ellipse plot 
    subplot(1,2,2); 
    plot(A_x, A_y); 
    axis equal; 
    xlabel('x'); 
    ylabel('y'); 
    title('ellipse plot');
    
    % add sigma_1 * u1, sigma_2 * u2 to ellipse plot 
    hold on; 
    plot([0, sigma_1 * u1(1)], [0, sigma_1 * u1(2)], 'LineWidth', 2);
    text(sigma_1 * u1(1), sigma_1 * u1(2), '  \sigma_1 u1', 'FontSize', 12);
    plot([0, sigma_2 * u2(1)], [0, sigma_2 * u2(2)], 'LineWidth', 2);
    text(sigma_2 * u2(1), sigma_2 * u2(2), '  \sigma_2 u2', 'FontSize', 12);
    hold off; 
end 

svd_plot(); 


% Q5


function [A,x,p, condA] = polyInterpOrApprox_mod(t, b, deg, nPlot, wantPlot, alg)
    % compare length of t and b
    if (length(t) ~= length(b))
        error("t and b are not the same length");
    end

    m = length(t); % if equal, set m to the length of t and b

    % compare m with deg + 1
    if (m < deg + 1)
        error("insufficient data, m < deg + 1")
    end

    % ensure valid choice for alg 
    if (alg <= 0 || alg > 6 || mod(alg,1) ~= 0)
        error("algorithm choice must be an integer between 1-6"); 
    end 

    % set up Vandermonde matrix with m rows, deg + 1 columns
    A = zeros(m, deg + 1);
    
    % we want the polynomial to be x(1) * t^deg + x(2) * t^(deg-1) + ... +
    % x(end) * t^(0), so we index j from 1 to deg + 1
    for j = 1 : deg + 1 % set the columns of A one at a time
        A(:,j) = t.^(deg + 1 - j); % componentwise exponentiation of vector t
    end
    
    % over-write variable condA in the case alg == 6
    condA = NaN; 

    switch alg
        case 1 
            x = A \ b;
        case 2
            x = (A'*A)\(A'*b);
        case 3
            % we take the transpose below to 
            % recover the lower triangular matrix G
            G = (chol(A'*A))'; 
            z = G \ A'*b; 
            x = G' \ z;
        case 4 
            [AA, q] = house(A); 
            [x, ~] = lsol(b, AA, q); 
        case 5
            [Q, R] = qr(A, 'econ'); 
            x = R \ Q'*b; 
        case 6
            [U,S,V] = svd(A, 'econ'); 
            x = V * (S \ U'*b); 
            % over-write condA
            condA = S(1,1)/S(end,end); 
    end 
    
    p = []; % initialize p as an empty list

    if (wantPlot ~= 0)
        % check if nPlot == 0
        if (nPlot == 0)
            figure; % clear current figure
        else 
            figure; 
            % create a grid with nPlot equally-spaced points from t(1) to t(n)
            tt = linspace(0, 1, nPlot);
    
            % compute the interpolating/approximating polynomial on the grid tt
            p = polyval(x, tt); 
            plot(tt, p); % plot polynomial p on tt
            hold on;
        end
    
        plot(t, b, 'ro') % plot original points as red circles
        
        % add axis labels and title
        xlabel('t') 
        ylabel('v')
        title(sprintf("Interpolating/Approximating Polynomial of deg %d using alg %d", deg, alg));
    end 
end


% helper functions : 

% A&G's house  
 function [A,q] = house(A)

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

% A&G's lsol  
function [x,nr] = lsol(b,AA,p)
    y = b(:); [m,n] = size(AA);
    % transform b
    for k=1:n
      u = [p(k);AA(k+1:m,k)];
      y(k:m) = y(k:m) - 2*u *(u'*y(k:m));
    end
    % form upper triangular R and solve
    R = triu(AA(1:n,:));
    x = R \ y(1:n); nr = norm(y(n+1:m));
end 

% test function 
function test_poly(t, b, deg, nPlot)

    format longG; 
    % allocate arrays to hold test results 
    algs = (1:6)'; 
    coeffs = cell(6, 1); 
    condAs = nan(6, 1); 

    for alg = 1:6

        [~, x, ~, condA] = polyInterpOrApprox_mod(t, b, deg, nPlot, 1, alg); 
        
        s = '';
        for i = 1:length(x)
            s = [s, sprintf('%22.16g', x(i))];  
        end
        coeffs{alg} = s; 

        if (~isnan(condA))
            condAs(alg) = condA; 
        end
    end 

    T = table(algs, coeffs, condAs, ...
        'VariableNames', {'algorithm', 'x', 'cond(A)'}); 

    writetable(T,'/Users/aryaman/numericalcomputinghw6.txt','Delimiter','\t','WriteRowNames',true);
    type /Users/aryaman/numericalcomputinghw6.txt
end

t = (linspace(0,1,101))'; 
b = cos(4*t); 
deg = 13; 
test_poly(t, b, deg, nPlot); 
