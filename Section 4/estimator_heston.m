function I = estimator_heston(dW,dt,Nq,TOL_Newton,option)
%Heston dynamics with FT-scheme

T = 1; r = 0.05; rho = -0.5; v0 = 0.04; theta = 0.0025; xi = 0.1; kappa = 1.0;

%Gaussian CDF
Phi = @(x)(1+erf(x/sqrt(2)))/2;

if option == 3 %Heston density
    X0 = 1; K = 1;
else
    X0 = 100; K = 100;
end

W = cumsum(dW,2);
dB = dW - dt/T*W(:,end);

%variance process
N = T/dt;
v = zeros(1,N-1);
vtemp = v0;

for j = 1:N-1
    vtemp = vtemp + kappa*(theta-max(vtemp,0))*dt + xi*sqrt(max(vtemp,0))*dW(1,j);
    v(j) = max(vtemp,0);
end
v = [v0 v];

%polynomial
G = @(x)  1 + r*dt + rho*sqrt(v).*dW(1,:) + sqrt(1-rho^2).*sqrt(v).*(dt/sqrt(T)*x + dB(2,:));
pol = @(x) X0*(prod(G(x),2)) - K;

diff =  sqrt(v).*dt/sqrt(T).*sqrt(1-rho^2);
dpol = @(x) X0*prod(G(x),2).*sum(diff./G(x),2);

a = Newtons_Method(pol,dpol,W(2,end),TOL_Newton);
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
    
elseif option == 3 % Heston density
    %Return estimate
    I = 1/sqrt(2*pi)*exp(-1/2*a.^2).*1./dpol(a);
    
elseif option == 4 % Heston Eurupean call Delta
    dXT = @(x) prod(G(x),2); 
    f = @(x) 1/2*(sign(pol(x))+1) .* 1/sqrt(2*pi).*exp(-x.^2/2) .* dXT(x);
    [x2,w2] = GaussLaguerre(Nq,alpha,a,b);
    %Return estimate
    I = exp(-r*T) .* sum(w2.*f(x2).*exp(x2-a));
    
elseif option == 5 % Heston European call Rho
    dW = dW(1,:);
    dB = dB(2,:); 
    f = @(x) exp(-r*T).*1/sqrt(2*pi).*exp(-x.^2/2).*(-T*max(pol(x),0) + 1/2*(sign(pol(x))+1).*dXRho(v,rho,dW,dB,X0,r,T,dt,N,x));
    [x1,w1] = GaussLaguerre(Nq,alpha,a,b);
    %Return estimate
    I = sum(w1.*f(x1).*exp(x1-a));  
    
elseif option == 6 % Heston digital call Delta
    dXT = @(x) prod(G(x),2);
    %Return estimate
    I = exp(-r*T).*1/sqrt(2*pi)*exp(-1/2*a^2)*dXT(a)*1/dpol(a); 
    
elseif option == 7 % Heston Digital call Rho
    %Return estimate
    I = exp(-r*T).*(-T + T*Phi(a) + 1/sqrt(2*pi)*exp(-1/2*a^2)*1/dpol(a).*dXRho(v,rho,dW(1,:),dB(2,:),X0,r,T,dt,N,a));
end

function ds = dXRho(v,rho,dW,dB,X0,r,T,dt,N,x)
F{1} = @(x) X0*dt;
F{2} = @(x) X0*(1 + r*dt + rho*sqrt(v(1)).*dW(1) + sqrt(1-rho^2).*sqrt(v(1)).*(dt/sqrt(T).*x + dB(1)));
J = F{2}(x);
ds = F{1}(x);
D = cell(1,N-1);
for k=2:N
    D{k} = @(x) 1 + r*dt + rho*sqrt(v(k)).*dW(k) + sqrt(1-rho^2).*sqrt(v(k)).*(dt/sqrt(T).*x + dB(k));
    ds = ds .* D{k}(x) + J.*dt;
    J = J .* D{k}(x);
end
end

end





