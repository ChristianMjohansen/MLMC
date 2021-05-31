function I = estimator_GBM(dW,dt,Nq,TOL_Newton,option)
%GBM dynamics with Euler-maruyama scheme

T = 1; r = 0.05; sig = 0.2;

if option == 3 %GBM Density
    X0 = 1; K = 1;
else
    X0 = 100; K = 100;
end

W = cumsum(dW);
dB = dW - dt/T*W(end);

G = @(x) 1 + r*dt + dt*sig/sqrt(T)*x + sig.*dB;
pol = @(x) X0*prod(G(x),2) - K;
dpol = @(x) X0*sig*dt/sqrt(T).*prod(G(x),2).*sum(1./G(x),2);

%Perform root finding
a = Newtons_Method(pol,dpol,W(end)/sqrt(T),TOL_Newton);

alpha = 0;
b = 1;
if option == 1 % European Call 
    f = @(x)  exp(-r*T) * max(pol(x),0) .* 1/sqrt(2*pi).*exp(-x.^2/2);
    % Pre-integration via Laguerre quadrature
    [x1,w1] = GaussLaguerre(Nq,alpha,a,b);
    %Return estimate
    I = sum(w1.*f(x1).*exp(x1-a));
    
elseif option == 2 % Digital option
    f = @(x) exp(-r*T) * 1/2*(sign(pol(x))+1) .* 1/sqrt(2*pi).*exp(-x.^2/2);
    % Pre-integration via Laguerre quadrature
    [x1,w1] = GaussLaguerre(Nq,alpha,a,b);
    %Return estimate
    I = sum(w1.*f(x1).*exp(x1-a));
    
elseif option == 3 % GBM Density
    %Return estimate
    I = 1/sqrt(2*pi)*exp(-1/2*a.^2).*1./dpol(a);
    
elseif option == 4 % GBM European call Delta
    dXT = @(x) prod(G(x),2);  
    f = @(x) 1/2*(sign(pol(x))+1) .* 1/sqrt(2*pi).*exp(-x.^2/2) .* dXT(x);
    [x2,w2] = GaussLaguerre(Nq,alpha,a,b);
    %Return estimate
    I = exp(-r*T) .* sum(w2.*f(x2).*exp(x2-a));
      
elseif option == 5 % GBM European call Vega  
    f = @(x) 1/2*(sign(pol(x))+1) .* 1/sqrt(2*pi).*exp(-x.^2/2).*dXsig(dB,X0,sig,r,T,dt,1/dt,x);
    [x1,w1] = GaussLaguerre(Nq,alpha,a,b);
    %Return estimate
    I = exp(-r*T) .* (sum(w1.*f(x1).*exp(x1-a)));

elseif option == 6 % GBM Digital call Delta
    Dn = @(x) 1 + r*dt + sig*(dt/sqrt(T)*x + dB);
    dXT = @(x) prod(Dn(x),2);
    I = exp(-r*T).*1/sqrt(2*pi)*exp(-1/2*a^2)*dXT(a)*1/dpol(a);     
    
elseif option == 7 % GBM Digital call Vega
    %Return estimate
    I = exp(-r*T).*1/sqrt(2*pi)*exp(-1/2*a^2)*1/dpol(a)*dXsig(dB,X0,sig,r,T,dt,1/dt,a);
end

function ds = dXsig(dB,X0,sig,r,T,dt,N,x)
F{1} = @(x) X0*(dt/sqrt(T).*x + dB(1));
F{2} = @(x) X0*(1 + r*dt + sig*(dt/sqrt(T).*x + dB(1)));
J = F{2}(x);
ds = F{1}(x);
D = cell(1,N-1);
for k=2:N
    D{k} = @(x) 1 + r*dt + sig*(dt/sqrt(T).*x + dB(k));
    ds = ds .* D{k}(x) + J.*(dt/sqrt(T).*x + dB(k));
    J = J .* D{k}(x);
end
end

end
