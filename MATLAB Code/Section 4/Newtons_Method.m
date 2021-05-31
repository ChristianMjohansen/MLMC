function x = Newtons_Method(f, df, x0, tol)
x  = x0;
%dx = 1e-5;
%df = @(x,dx) (f(x + dx) - f(x - dx))/(2*dx);
while abs(f(x)) > tol
    x = x - f(x)./df(x);
    %x = x - f(x)./df(x,dx);
end
end
