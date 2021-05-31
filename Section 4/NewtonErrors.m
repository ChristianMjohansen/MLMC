%Newtons method errors
clc; clear;  T = 1; N = 2^6; dt = T/N; sig = 0.2; X0 = 100; K = 100; r = 0.05; 

%increments
dW = bb(randn(N,1),T)';
W = cumsum(dW,2);
dB = dW-dt/T.*W(end);

%polynomial
F = @(x) 1 + r*dt + dt*sig/sqrt(T).*x + sig.*dB;
pol = @(x) X0*prod(F(x),2) - K;
dpol = @(x) X0*sig*dt/sqrt(T)*prod(F(x)).*sum(1./F(x),2);

%Newton's method
x = W(end)/sqrt(T);
TOL = 1e-12;

[a,n,err] = Newtons_Method2(pol,dpol,x,TOL);

semilogy(1:n,err,'r-*');
grid on;

xlabel('Number of iterations')
ylabel('Absolute error')
axis([1 length(err) min(err) max(err) ])


print(gcf,'-depsc','-painters','Results/NewtonErrors.eps')
