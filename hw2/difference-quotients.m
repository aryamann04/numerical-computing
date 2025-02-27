function differenceQuotientErrors(print_details)
%
% translate Program 11.1 from Ch. 11 of the book to Matlab with 
% more detailed output, showing how the cancellation occurs,
% and generate the log-log plot discussed on p.77
%
% for the function f(x)=sin(x), for which the derivative is f'(x)=cos(x),
% compute the difference quotient (f(x+h) - f(x))/h for various h and
% compare this with the true derivative. The difference is the error.
% The error is dominated by DISCRETIZATION ERROR (OR TRUNCATION ERROR)
% when h is large and by CANCELLATION ERROR (ROUNDING ERROR) when h is small. 
% Use x = 1.
% 
% default input is 1
%
if nargin == 0
    print_details = 1;
end
%
% define the anonymous function f and its true derivative
%
f = @(x)sin(x);
fprime = @(x)cos(x);
x = 1;
fx = f(x);
derivx = fprime(x);
% print f(x) to full precision since it cancels with f(x+h) when h is too small
fprintf('Using formulas,   f(x) = %20.16e   and   f''(x) = %10.7e\n',fx,derivx)
% collect all the values of h in a vector h using the .^ operator
h_all = 10.^(-20:0); % vector with components 10^-20, 10^-19, ..., 1, 0
if print_details % print column headers
    fprintf('  h            x+h             f(x+h)             f(x+h)-f(x)  (f(x+h)-f(x))/h  error\n')
end
for k=1:length(h_all)
    h = h_all(k);
    xph = x + h;
    fxph = f(xph);
    dif = fxph - fx;
    dif_quo = dif/h;
    error = abs(derivx - dif_quo);
    if print_details % print fxph to full precision because it is cancelling with fx
        fprintf('%6.2e  %11.7e  %20.16e  %11.7e  %11.7e  %6.2e\n',...
            h, xph, fxph, dif, dif_quo, error)
    end
    error_all(k) = error;
end
loglog(h_all,error_all,'*')
xlabel('h')
ylabel('error');
title('errors in approx of derivative of sin at x=1 by difference quotients')

