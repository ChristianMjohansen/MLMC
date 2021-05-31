function V = digital_call(r,sigma,T,X,K,opt)  

d1 = ( log(X) - log(K) + (r+1/2*sigma^2)*T )/(sigma*sqrt(T));
d2 = ( log(X) - log(K) + (r-1/2*sigma^2)*T )/(sigma*sqrt(T));

%Gaussian CDF
Phi = @(x)(1+erf(x/sqrt(2)))/2;

switch opt
  case 'value'
    V = exp(-r*T).*Phi(d2);

  case 'delta'
    V = exp(-r*T)*(1/(sqrt(2*pi)).*exp(-1/2*d2.^2))./(sigma*sqrt(T)*X);
    
  case 'vega'
    V = -exp(-r*T)*(d1./(sigma*sqrt(2*pi)).*exp(-1/2*d2.^2));

end



