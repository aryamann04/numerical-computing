function w = baryConstruct(x)
    % force x into a column vector 
    x = x(:); 
    n = numel(x); 
    w = ones(n, 1); 
    
    % shift indices up by 1
    for j = 1 : n
        % w(j) = 1 / (vectorized product skipping x(j))
        w(j) = 1 ./ ( prod(x(j) - x([1:j-1, j+1:n])) );
    end 
end

function p = baryEval(w, x, y, z)
    
    sz = size(z); 

    % force to column vectors for consistency 
    w = w(:); x = x(:); y = y(:); z = z(:); 

    % initialize vector to store the evaluations p(z_j)
    p = zeros(size(z)); 
    
    % tolerance for equality check in loop 
    tol = eps(max(abs(x))); 

    % evaluate at each z_j  
    for j = 1 : numel(z)
        
        % check if z_j == x_i for any i
        % note that to avoid dividing by an extremely small number, we set 
        % p(z_j) = y(i) if |x(i) - z(j)| < tol = ulp(max(|x_i|))
        
        % we have to first find a matching index i s.t. z(j) = x(i) 
        i = find(abs(z(j) - x) < tol, 1); 
        if (~isempty(i))
            % z(j) = x(i) so p(z(j)) = y(i)
            p(j) = y(i); 
        else 
            d = 1 ./ (z(j) - x); % 1/(z(j) - x) vectorized

            % numerator = sum k=1:n w(k)y(k)/(z(j)-x(k))
            numerator = sum( w .* y .* d ); 

            % denominator = sum k=1:n w(k)/(z(j)-x(k))
            denominator = sum (w .* d); 
            
            p(j) = numerator / denominator; 
        end 
    end

    p = reshape(p, sz); 
    z = reshape(z, sz); 
end

% Runge function, accepts a column vector x as an argument 
function y = f(x)

    y = 1 ./ (1 + 25 * x.^2);

end 

function q1_main()

    n = 10; % 11 points --> degree 10 polynomial 

    % x: column vector of 11 equally spaced points in [-1,1]
    x = linspace(-1, 1, n + 1).'; 

    % y = f(x), where f is the Runge function 
    y = f(x); 

    % z: column vector of N + 1 equally spaced points, N >> n
    N = 10 * n; % set to 10 * n for simplicity 
    z = linspace(-1, 1, N + 1).'; 

    % call baryConstruct and baryEval
    w = baryConstruct(x); 
    p = baryEval(w, x, y, z); 

    % compare to MATLAB's natural cubic spline 
    sz = spline(x, y, z); 

    % plot for part (b)
    figure(1); hold on; grid on;
    
    % (x(i), y(i))
    plot(x, y, 'r*', 'DisplayName', 'interpolation points (x(i), y(i))');

    % exact Runge function
    plot(z, f(z), 'k-', 'DisplayName', 'Runge f(x)');
    
    % barycentric interpolant
    plot(z, p, 'b-', 'DisplayName', 'barycentric (degree 10)');

    % MATLAB's cubic spline 
    plot(z, sz, 'g-', 'DisplayName', 'cubic spline');

    legend('Location','SouthEast');
    title('Runge function: barycentric vs spline on [-1,1]');
    xlabel('x'); ylabel('y');
    hold off;

    % part (c) : plotting the individual cubics   

    pp = spline(x, y); 

    % 10 different RGB values for plot colors 
    colors = lines(10); 
    
    figure(2); hold on; grid on;
    for i = 1 : size(pp.coefs, 1)
        % get x_i, x_{i+1}
        x_i = pp.breaks(i); 
        x_ip1 = pp.breaks(i+1); 

        % boolean vector so we can color the relevant segment 
        % for z in [x_i, x_{i+1}] black
        in_interval = (z >= x_i) & (z <= x_ip1); 
        
        p_i = polyval(pp.coefs(i, :), z - x_i);

        plot(z, p_i, '-', 'Color', colors(i, :), ...
            'DisplayName', sprintf('cubic #%d', i)); 

        % overwrite part of plot that is part of the spline 
        plot(z(in_interval), p_i(in_interval), '--k', 'LineWidth', 1.2, ...
            'HandleVisibility', 'off'); 
    end 

    % add the exact Runge function 
    plot(z, f(z), 'k-', 'DisplayName', 'Runge f(x)');
    title('Runge function: piecewise cubics (solid) and spline segments (dashed)');
    xlabel('x');  ylabel('y');
    legend('Location','SouthEastOutside');
    hold off;
end

% q1_main(); 

function x = newton(f, fderiv, x0, nsteps)
    
    % print out initial guess 
    fprintf("iteration 0:\tx0 = %23.16e\tf(x0) = %g\n", x0, f(x0)); 
    
    x = x0; 
    for k = 1 : nsteps
        x = x - ( f(x) / fderiv(x) ); 
        fprintf("iteration %d:\tx%d = %23.16e\tf(x%d) = %g\n", ...
            k, k, x, k, f(x));
    end
    fprintf("\n\n"); 
end

function q2_main()

    % test 1: f1(x) = log(x) - sin(x)
    f1 = @(x) log(x) - sin(x); 
    f1deriv = @(x) 1./x - cos(x); 
    x1 = newton(f1, f1deriv, 1, 10); 

    % test 2: f2(x) = x^3 
    f2 = @(x) x.^3; 
    f2deriv = @(x) 3*x.^2; 
    x2 = newton(f2, f2deriv, 0.1, 10); 

    % test 3: f3(x) = sin(x)
    f3 = @(x) sin(x); 
    f3deriv = @(x) cos(x); 
    x3 = newton(f3, f3deriv, 3, 10); 

end

% q2_main(); 

function [v, convg, r] = boundary_ode(alpha, lambda, plot_on)
    
    % default to lambda = 1, plotting off 
    if nargin < 2, lambda = 1; end
    if nargin < 3, plot_on = false; end
    
    N = 100; 
    h = 1/N; 
    max_iter = 50; 
    tol = 10^(-10); 
    convg = false; 
    
    % interior t-values (t_1,...,t_{n-1})
    t = (1:N-1).' * h; 
    
    % initial guess: v0 = alpha * (t1(1-t1), ... , tn(1-tn))'
    v = alpha * t .* (1 - t); 

    % initialize vector iterates to track v for plotting
    iterates = v';
    
    % print table header
    fprintf("-----------------------------\n");
    fprintf("%3s    %20s\n", "k", "||g(v^{(k)})||_2"); 
    fprintf("-----------------------------\n");
    for k = 0 : max_iter

        % discretize the ODE by Taylor's approximation for v''(t)

        % v_{i-1} : right-shift v and set v0 = 0
        vim1 = ([0; v(1 : end-1)]); 
        % v_{i+1} : left-shift v and set vN = 0
        vip1 = ([v(2 : end); 0]); 
        % construct the column vector g(v)
        g = (vip1 - 2*v + vim1)/(h^2) + lambda * exp(v); 

        % check for convergence and print residual 
        r = norm(g, 2); 
        fprintf('%3d    %20.12e\n', k, r); 

        if (r < tol)
            convg = true; 
            break; 
        end 

        % set up Jacobian matrix for the discretized ODE: 
        %   super diagonal  : 1/h^2             = (1/h^2) * (1)
        %   main diagonal   : -2/h^2 + e^{v(i)} = (1/h^2) * (-2 + h^2e^{v(i)})
        %   sub diagonal    : 1/h^2             = (1/h^2) * (1)
        
        % note: the only difference when lambda =/= 1 is the constant factor 
        %       multiplying exp(v) in the main diagonal 

        main = -2/(h^2) + lambda * exp(v); 
        off = (1/h^2) * ones(N-2, 1); 
        J = spdiags([[0; off], main, [off; 0]], -1:1, N-1, N-1); 

        % solve for the direction vector p and execute newton's step
        p = J \ -g; 
        v = v + p; 

        % track current iterate
        iterates(end + 1, :) = v';
    end

    fprintf("-----------------------------\n");

    % check if algorithm converged 
    if ~convg 
        fprintf("algorithm did not converge: norm(g) = %23.16e\n", r);
    else
        fprintf("converged in %d iterations\n", k);
    end 
    
    % return v with endpoints included 
    v = [0; v; 0];

    % plot solution and iterates
    if plot_on
        % expand t to also include endpoints 
        t = [0; t; 1]; 
        figure; hold on; grid on; 

        % get darker as iterates converge to solution 
        colors = [linspace(0.8, 0, k)', ... 
            linspace(0.8, 0, k)', ... 
            linspace(0.8, 0, k)'];  
    
        % plot iterates
        for i = 1 : k 
            plot(t, [0; iterates(i, :)'; 0], 'Color', colors(i, :), ...
                'LineWidth', 1.0, ...
                'DisplayName', sprintf('iter %d', i));
        end

        % plot solution v(t)
        plot(t, v, '-k', 'LineWidth', 2, 'DisplayName', 'sol v(t)');
        
        title(sprintf("newton iterates and solution v(t) with " + ...
            "\lambda = %g and \alpha = %g", lambda, alpha));
    
        xlabel('t'); ylabel('v(t)'); 
        legend show; hold off; 
    end
end

function find_critical_lambda(lambdas, alphas)
    
    results = [];            
    lambda_vals = [];        
    v_mid_vals = [];

    for lambda = lambdas
        v_midpoints = [];    % stores distinct v(0.5) values

        for alpha = alphas
            try
                [v, convg, ~] = boundary_ode(alpha, lambda, false);

                if convg

                    v_mid = v(round(end/2));
                    
                    lambda_vals(end+1) = lambda;
                    v_mid_vals(end+1) = v_mid;

                    % check if v_mid is distinct
                    if all(abs(v_mid - v_midpoints) > 1e-2)
                        v_midpoints(end+1) = v_mid;
                    end
                end
                
            catch
                % ignore errors 
            end
        end

        results(end+1,:) = [lambda, length(v_midpoints)];
    end

    figure;
    subplot(1,2,1);
    plot(results(:,1), results(:,2), '-o', 'LineWidth', 1.5);
    xlabel('\lambda'); ylabel('number of solutions');
    title('count of distinct solutions');
    grid on;

    subplot(1,2,2);
    scatter(lambda_vals, v_mid_vals, 30, 'filled');
    xlabel('\lambda'); ylabel('v(0.5)');
    title('\lambda vs v(0.5)');
    grid on;
end

function q3_main()

    % part (b) 
    
    % solution 1: alpha = 0 
    boundary_ode(0, 1, true);

    % solution 2: alpha = 20
    boundary_ode(20, 1, true); 
    
    % part (c) 
    
    % experimenting with increasing values of lambda 

    lambdas = 1:0.5:5;
    alphas = [1, 20];

    figure; hold on; grid on;
    colors = lines(length(lambdas));
    legend_handles = gobjects(length(lambdas), 1);
    legend_labels = strings(length(lambdas), 1);

    for i = 1:length(lambdas)
        lambda = lambdas(i);
        color_i = colors(i, :);

        plotted = false;

        for alpha = alphas

            [v, convg, ~] = boundary_ode(alpha, lambda, false);

            if convg
                t = linspace(0, 1, length(v));
                h = plot(t, v, 'Color', color_i, 'LineWidth', 1.5);

                if ~plotted
                    legend_handles(i) = h;
                    legend_labels(i) = sprintf('\\lambda = %.1f', lambda);
                    plotted = true;
                end
            else
                fprintf("did not converge for \lambda = %.2f " + ...
                    "with \alpha = %.2f\n", lambda, alpha);
            end
        end
    end

    title('solutions v(t) for various \lambda');
    xlabel('t'); ylabel('v(t)');

    legend(legend_handles(isgraphics(legend_handles)), ...
        legend_labels(isgraphics(legend_handles)), 'Location', 'best');

    hold off;
    
    % finding the critical lambda 

    lambdas2 = 1 : 0.1 : 6;
    alphas2 = 1 : 2 : 50;
    
    find_critical_lambda(lambdas2, alphas2); 
    
    % result: 3.4 < lambda < 3.6

    lambdas3 = 3.4 : 0.001 : 3.6;

    find_critical_lambda(lambdas3, alphas2); 

    % result: 3.503 < lambda < 3.520

end 

% q3_main(); 

boundary_ode(0, 1, true);
boundary_ode(20, 1, true);