function [A,x,p] = polyInterpOrApprox(t, b, deg, nPlot)
    % compare length of t and b
    if (length(t) ~= length(b))
        error("t and b are not the same length");
    end

    m = length(t); % if equal, set m to the length of t and b

    % compare m with deg + 1
    if (m < deg + 1)
        error("insufficient data, m < deg + 1")
    end 

    % set up Vandermonde matrix with m rows, deg + 1 columns
    A = zeros(m, deg + 1);
    
    % we want the polynomial to be x(1) * t^deg + x(2) * t^(deg-1) + ... +
    % x(end) * t^(0), so we index j from 1 to deg + 1
    for j = 1 : deg + 1 % set the columns of A one at a time
        A(:,j) = t.^(deg + 1 - j); % componentwise exponentiation of vector t
    end
    
    x = A \ b; % Matlab's backslash operator solves Ax=b
    
    % check if nPlot == 0
    if (nPlot == 0)
        p = []; % initialize p as an empty list
        clf; % clear current figure
    else 
        % create a grid with nPlot equally-spaced points from t(1) to t(n)
        tt = linspace(0, 20, nPlot);

        % compute the interpolating/approximating polynomial on the grid tt
        p = polyval(x, tt); 

        clf % clear current figure
        plot(tt, p); % plot polynomial p on tt
        hold on
    end

    plot(t, b, 'ro') % plot original points as red circles
    
    % add axis labels and title
    xlabel('t') 
    ylabel('v')
    title('interpolating/approximating polynomial')
end