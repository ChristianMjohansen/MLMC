function out = basket_call(r,sigma,T,X0,K,rho,d,opt)

M = 10^7; % number of samples
O = eye(d) + rho*(ones(d)-eye(d)); %chol decomposition
L = chol(O,'lower');
W = sqrt(T)*L*randn(d,M);
X = X0.*exp((r-0.5*sigma^2)*T + sigma*W);
X = 1/d*sum(X); % average
switch opt
    case 'Euro call'
        F1 = exp(-r*T)*max(X-K,0); % call option
        out = sum(F1)/M;
    case 'digi call'
        F2 = exp(-r*T)*1/2*(sign(X-K)+1); %digital call
        out = sum(F2)/M;
end


