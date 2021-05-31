function out = heston(r,kappa,theta,xi,rho,X0,K,v0,T,opt)
%Routine for estimating European, digital options, risk neutral density, and a few Greeks of the Heston model.

int1 = @(x,X0,v0,theta,kappa,xi,r,rho,T,K)... 
real(exp(-1i.*x*log(K)).*charFkt(x-1i,r,kappa,theta,xi,rho,X0,v0,T)./(1i*x.*charFkt(-1i,r,kappa,theta,xi,rho,X0,v0,T))); 
int1 = integral(@(x)int1(x,X0,v0,theta,kappa,xi,r,rho,T,K),0,Inf); % numerical integration
P1 = 1/2+int1/pi; 

int2 = @(x,X0,v0,theta,kappa,xi,r,rho,T,K) real(exp(-1i.*x*log(K)).*charFkt(x,r,kappa,theta,xi,rho,X0,v0,T)./(1i*x));
int2 = integral(@(x)int2(x,X0,v0,theta,kappa,xi,r,rho,T,K),0,Inf); % numerical integration
int2 = real(int2);
P2 = 1/2+int2/pi;

int3 = @(x,X0,v0,theta,kappa,xi,r,rho,T,K) real(exp(-1i.*x*log(K)).*charFkt(x,r,kappa,theta,xi,rho,X0,v0,T));
int3 = integral(@(x)int3(x,X0,v0,theta,kappa,xi,r,rho,T,K),0,Inf); % numerical integration
int3 = real(int3);
dPX = int3/(pi*X0);

int4 = @(x,X0,v0,theta,kappa,xi,r,rho,T,K) real(exp(-1i.*x*log(K)).*T.*charFkt(x,r,kappa,theta,xi,rho,X0,v0,T));
int4 = integral(@(x)int4(x,X0,v0,theta,kappa,xi,r,rho,T,K),0,Inf); % numerical integration
int4 = real(int4);
dPr = int4/pi;

%
% Integrals for estimating Heston density, using Breeden and Litzenberger
% formula.
% 
dint1 = @(x,X0,v0,theta,kappa,xi,r,rho,T,K)...
real((((1i.*x-x.^2).*K.^(-1i.*x-2)).*charFkt(x-1i,r,kappa,theta,xi,rho,X0,v0,T))./(1i*x.*charFkt(-1i,r,kappa,theta,xi,rho,X0,v0,T))); 
dint1 = integral(@(x)dint1(x,X0,v0,theta,kappa,xi,r,rho,T,K),0,Inf); % numerical integration
p1 = dint1/pi; 

dint2 = @(x,X0,v0,theta,kappa,xi,r,rho,T,K) real((-K.^(-1i.*x-1)).*charFkt(x,r,kappa,theta,xi,rho,X0,v0,T));
dint2 = integral(@(x)dint2(x,X0,v0,theta,kappa,xi,r,rho,T,K),0,Inf); % numerical integration
dint2 = real(dint2);
p2 = dint2/pi; 

dint3 = @(x,X0,v0,theta,kappa,xi,r,rho,T,K) real(((1i.*x-x.^2) .*K.^(-1i.*x-2)).*charFkt(x,r,kappa,theta,xi,rho,X0,v0,T)./(1i*x));
dint3 = integral(@(x)dint3(x,X0,v0,theta,kappa,xi,r,rho,T,K),0,Inf); % numerical integration
dint3 = real(dint3);
p3 = dint3/pi; 

% European, digital call options, various greeks and density.
switch opt
    case 'Euro call'
        out = X0*P1-exp(-r*T)*K*P2; %put = call-X0+K*exp(-r*T); 
    case 'digi call'
        out = exp(-r*T)*P2;
    case 'Euro delta'
        out = P1;
    case 'Euro rho'
        out = K*exp(-r*T)*T*P2;
    case 'digi delta'
        out = exp(-r*T)*dPX;
    case 'digi rho'
        out = -T*exp(-r*T)*P2 + exp(-r*T)*dPr;
    case 'density'
        out = X0*exp(r*T)*p1 - 2*p2 - K*p3;
end

%Characteristic function
function out = charFkt(x,r,kappa,theta,xi,rho,X0,v0,T)
a = -x.*x/2 - 1i*x/2;
b = kappa-rho*xi*1i*x;
c = xi^2/2;
M = sqrt(b.*b-4*a*c);
zplus = (b+M)/xi^2;
zminus = (b-M)/xi^2;
N = zminus./zplus;
A = kappa*(zminus*T-(2/xi^2).*log((1-N.*exp(-M*T))./(1-N)));
D = zminus.*(1-exp(-M*T))./(1-N.*exp(-M*T));
out = exp(A*theta+D*v0+1i*x*log(X0*exp(r*T))); 




