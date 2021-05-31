function V = european_call(r,sigma,T,X,K,opt)

d1 = (log(X) - log(K) + (r+1/2*sigma^2)*T )/(sigma*sqrt(T));
d2 = (log(X) - log(K) + (r-1/2*sigma^2)*T )/(sigma*sqrt(T));

%Gaussian CDF
Phi = @(x)(1+erf(x/sqrt(2)))/2;

switch opt
  case 'value'
    V = X.*Phi(d1) - exp(-r*T)*K.*Phi(d2);

  case 'delta'
    V = Phi(d1);

  case 'vega'
    V = X.*sqrt(T)./sqrt(2*pi)*exp(-1/2*d1.^2);
 
end
