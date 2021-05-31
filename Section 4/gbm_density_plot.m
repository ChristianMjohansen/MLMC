%Monte Carlo with numerical smoothing (GBM density)
clc; clear;
T = 1; X0 = 1; r = 0.05; sig = 0.2;  N = 2^8; dt = T/N; M = 5000; 

TOL_Newton = 1e-3;
npoints = 50;
K = linspace(0.3,2,npoints);
res = zeros(1,M);
Y = zeros(1,length(K));

for i=1:length(K)
    parfor j=1:M
        dW = sqrt(dt)*randn(1,N);
        W = cumsum(dW,2);
        dB = dW - dt/T*W(end);
        
        %polynomial
        G = @(x) 1 + r*dt + dt*sig/sqrt(T)*x + sig.*dB;
        pol = @(x) X0*prod(G(x),2) - K(i);

        dpol = @(x) X0*sig*dt/sqrt(T)*prod(G(x),2).*sum(1./G(x));
        
        %Perform root finding
        a = Newtons_Method(pol,dpol,W(end)/sqrt(T),TOL_Newton);
        
        res(j) = exp(-1/2*a.^2)*1/dpol(a);
    end
    Y(i) = 1/sqrt(2*pi)*mean(res);
end

%reference
den = @(u) 1./(sqrt(2*pi*T).*sig*u).*exp(-1/2*(log(u)-log(X0)-(r-sig^2/2)*T).^2./(sig^2*T));

plot(K,Y,'b*',K,den(K),'--+r') 
xlabel('K'); ylabel('\rho_X');
legend('GBM density (smoothed)','Reference','Fontsize',10)
print(gcf,'-depsc','-painters',['Results/GBM_den.eps'])


