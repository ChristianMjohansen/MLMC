function I = estimator_4D_GBM(dW,dt,A,Ainv,L,Nq,TOL_Newton,option)
%4D GBM dynamics with Euler-maruyama scheme

T = 1; sig = 0.2; X0 = 100; K = 100; r = 0.05; 

W = cumsum(dW,2);
dB = dW - dt/T*W(:,end);
%linear mapping
Y = A*W(:,end)/sqrt(T);

%multivariate GBM
dif = @(j,x) sig*(dt/sqrt(T)*(Ainv(j,1)*x + sum(Ainv(j,2:end)*Y(2:end))) + dB(j,:));
drift = 1 + r*dt;
F1 = @(x) drift + L(1,1)*dif(1,x);
F2 = @(x) drift + L(2,1)*dif(1,x) + L(2,2)*dif(2,x);
F3 = @(x) drift + L(3,1)*dif(1,x) + L(3,2)*dif(2,x) + L(3,3)*dif(3,x);
F4 = @(x) drift + L(4,1)*dif(1,x) + L(4,2)*dif(2,x) + L(4,3)*dif(3,x) + L(4,4)*dif(4,x);

%polynomial
pol = @(x) 1/4*X0.*(prod(F1(x),2) + prod(F2(x),2) + prod(F3(x),2) + prod(F4(x),2)) - K;

dpol = @(x) 1/4*X0*sig*dt/sqrt(T)*(prod(F1(x),2)*sum(1./F1(x))*Ainv(1,1) + prod(F2(x),2)*sum(1./F2(x))*Ainv(2,1)...
                 + prod(F3(x),2)*sum(1./F3(x))*Ainv(3,1) + prod(F4(x),2)*sum(1./F4(x))*Ainv(4,1));

%Perform root finding
a = Newtons_Method(pol,dpol,Y(1),TOL_Newton);

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
end

end