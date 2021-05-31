function [x,n,err] = Newtons_Method(f, df, x0, tol)
x  = x0;
n = 0;
%dx = 1e-6;
%df = @(x,dx) (f(x + dx) - f(x - dx))/(2*dx);
err = [];
while abs(f(x)) > tol
    err = [ err  abs((x - f(x)/df(x))-x) ];
    x = x - f(x)/df(x);
    n = n + 1;
end
end
