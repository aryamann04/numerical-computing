function plotPieceWiseLinear(f, x, z)

    n = length(x); 
    figure; grid on; hold on; 
    
    for i = 1: n-1
        
        interval = (z >= x(i)) & (z <= x(i+1));
        y = ((f(x(i+1)) - f(x(i)))/(x(i+1) - x(i))) * (z(interval) - x(i)) + f(x(i)); 
        plot(z(interval), y, 'k-', 'DisplayName', 'interpolant'); 

    end 

    y = arrayfun(f, z); 
    plot(z, y, 'r-', 'DisplayName', 'f'); 
    xlabel('x'); ylabel('y'); 
    title('linear interploation of f(x)'); 
    legend show; hold off; 

end

function [afinal, bfinal] = bisection(f, a, b, tol)

    if a > b 
        error("a > b"); 
    elseif tol <= 0 
        error("tol <= 0"); 
    elseif f(a)*f(b) > 0 
        error("f(a) * f(b) > 0"); 
    end 

    while(true)
        
       if b - a <= tol
           break;
       end 

       mid = (a+b)/2; 

       if f(mid) * f(b) <= 0
           a = mid; 
       else 
           b = mid; 
       end 
    end 

    afinal = a; 
    bfinal = b; 

end 

function [mu, v] = geteig(A, tol_big, tol_small, maxit) 

    if norm(A - A') > 0
        error("A =\= A^T"); 
    elseif tol_big <= tol_small 
        error("tol_big <= tol_small");
    elseif tol_small <= 1e-15
        error("tol_small <= 1e-15"); 
    end 

    [n, ~] = size(A); 
    v = ones(n,1); v = v / norm(v);
    rayleigh = false; 

    for k = 1 : maxit 

        if rayleigh == true
            v_hat = (A - eye(n) * mu) \ v; 
        else 
            v_hat = A * v; 
        end 

        v = v_hat / norm(v_hat); 
        mu = v' * A * v; 
        r = norm(A * v - mu * v); 
        fprintf("iter %d, resid = %20.16e\n", k, r); 

        if r <= tol_big && rayleigh == false 
            rayleigh = true; 
            fprintf("\n****switching to rayleigh iteration (iter %d onward)****\n\n", k+1);
        end 
        
        if r <= tol_small, break; end 

    end 
    
    fprintf("\neigenvalue = %15.9e, eigenvector = \n", mu);
    disp(v); 
end 

function test()

    % Q1

    x = linspace(0,10, 11); 
    z = linspace(0, 10, 1001); 
    f = @(x) (x-10)*(x-2); 

    % plotPieceWiseLinear(f, x, z); 

    % Q2

    a = 1; b = 3; tol = 10^(-6);
    % [a, b] = bisection(f, a, b, tol); 
    % fprintf("zero in [%20.16e, %20.16e]\n", a, b); 

    % Q3
    A = [2 1 0;
     1 3 1;
     0 1 2];
    geteig(A, 1e-2, 1e-8, 100);

end

test(); 