%Gauss lageurre errors
clc; clear; format long;
alpha = 0;
b = 1;
Nq = 1;
pole = 0;

%Integrand 
f = @(x) 1/2*(sign(x)+1).*exp(-1/2*x.^2);
%f = @(x) max(x,0)*1/sqrt(2*pi).*exp(-1/2*x.^2);

%right side integral
[x1,w1] = GaussLaguerre(Nq,alpha,pole,b);
res = sum(w1.*f(x1).*exp(x1-pole));

TOL = 1e-12;
eps = 1; 
err=[]; 
while eps > TOL
    Nq = Nq + 1;
    [x1,w1] = GaussLaguerre(Nq,alpha,pole,b);
    y = sum(w1.*f(x1).*exp(x1-pole));
    eps = abs(y - res);
    err = [err eps];
    res = y;
end

semilogy(1:length(err),err)
grid on
axis([0 length(err) 0 1])
xlabel('N_q')
ylabel('Absolute error')

print(gcf,'-depsc','-painters','Results/LaguerreErrors.eps')














