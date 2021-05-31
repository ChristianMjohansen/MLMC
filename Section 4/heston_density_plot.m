%Monte Carlo with numerical smoothing (Heston density)
clc; clear;
T = 1; X0 = 1; r = 0.05; rho = -0.5; v0 = 0.04; theta = 0.0025; xi = 0.1; kappa = 1.0;

N = 2^8; dt = T/N; TOL_Newton = 1e-3; M = 50000;

res = zeros(1,M);
npoints = 50;
K = linspace(0.4,1.7,npoints);
temp1 = zeros(1,length(K)); 
temp2 = zeros(1,length(K));

for i=1:length(K)
    parfor j=1:M
        dW = sqrt(dt)*randn(2,N);
        W = cumsum(dW,2);
        dB = dW - dt/T*W(:,end);
        
        %Heston with FT-scheme
        v = zeros(1,N-1);
        vtemp = v0;
        for k = 1:N-1
            vtemp = vtemp + kappa*(theta-max(vtemp,0))*dt + xi*sqrt(max(vtemp,0))*dW(1,k);
            v(k) = max(vtemp,0);
        end
        v = [v0 v];
   
        %polynomial
        G = @(x)  1 + r*dt + rho*sqrt(v).*dW(1,:) + sqrt(1-rho^2).*sqrt(v).*(dt/sqrt(T)*x + dB(2,:));
        pol = @(x) X0*prod(G(x),2) - K(i);
        
        diff =  sqrt(v).*sqrt(1-rho^2).*dt/sqrt(T);
        dpol = @(x) X0*prod(G(x),2).*sum(diff./G(x));
        
        %Perform root finding
        a = Newtons_Method(pol,dpol,W(2,end)/sqrt(T),TOL_Newton);

        res(j) = exp(-1/2*a.^2).*1./dpol(a);
    end
    %smoothed density
    temp1(i) = 1/sqrt(2*pi)*mean(res);
    
    %reference density
    temp2(i)= heston(r,kappa,theta,xi,rho,X0,K(i),v0,T,'density');
end
plot(K,temp1,'b*',K,temp2,'--+r')
xlabel('K'); ylabel('\rho_X');
legend('Heston density (smoothed)','Reference','Fontsize',10)
print(gcf,'-depsc','-painters',['Results/Heston_den.eps'])



