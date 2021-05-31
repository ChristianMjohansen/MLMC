%Level estimator for standard MLMC
function [ sums, cost] = level_estimator(l,M,model,option)

T = 1; X0 = 100; K = 100; r = 0.05; sig = 0.2; rho = -0.5; v0 = 0.04; theta = 0.0025; xi = 0.1; kappa = 1.0;

sums(1:6) = 0;

nf = 2^l;
nc = nf/2;
hf = T/nf;
hc = T/nc;

Xf = X0*ones(1,M);
vf = v0*ones(1,M); 
Xc = Xf;
vc = vf;

dX0f = ones(1,M);
dX0c = dX0f;

dSigf = zeros(1,M);
dSigc = dSigf;

dRhof= zeros(1,M);
dRhoc = dRhof;

 
%Gaussian CDF
Phi = @(x)(1+erf(x/sqrt(2)))/2;

if model == 1  %GBM, euler-maruyama scheme
    alf = 1;
    if l==0 % Just one timestep
        dWf = sqrt(hf)*randn(1,M);
        Df = 1 + r*hf + sig.*dWf;
        dSigf = dSigf.*Df + Xf.*dWf;
        dX0f = dX0f.*Df;
        Xf = Xf.*Df;
    else
        for k = 1:nc
            dWf = sqrt(hf)*randn(2,M);
            for j = 1:2
                Df = 1 + r*hf + sig.*dWf(j,:);
                dSigf = dSigf.*Df + Xf.*dWf(j,:);
                dX0f = dX0f.*Df;
                Xf = Xf.*Df;
            end
            dWc  = dWf(1,:) + dWf(2,:);% Construct coarse brownian increments from the fine ones
            Dc = 1 + r*hc + sig.*dWc;
            dSigc = dSigc.*Dc + Xc.*dWc;
            dX0c = dX0c.*Dc;
            Xc = Xc.*Dc;
        end
    end
     
elseif model == 2 %Heston full truncation scheme
    alf = 1;
    if l==0  % Just one timestep
        dWf1 = sqrt(hf)*randn(1,M);
        dWf2 = rho*dWf1+sqrt(1-rho^2)*sqrt(hf)*randn(1,M);
        Df = 1 + r*hf + sqrt(max(vf,0)).*dWf2;
        dRhof = dSigf.*Df + Xf.*hf;
        dX0f = dX0f.*Df;
        Xf = Xf.*Df;
        vf = vf + kappa*(theta-max(vf,0))*hf + xi*sqrt(max(vf,0)).*dWf1;
    else
        for k = 1:nc
            dWf1 = sqrt(hf)*randn(2,M);
            dWf2 = rho*dWf1 + sqrt(hf)*randn(2,M)*sqrt(1-rho^2);
            for j = 1:2
                Df = 1 + r*hf + sqrt(max(vf,0)).*dWf2(j,:);
                dRhof = dRhof.*Df + Xf.*hf;
                dX0f = dX0f.*Df;
                Xf = Xf.*Df;
                vf = vf + kappa*(theta-max(vf,0))*hf + xi*sqrt(max(vf,0)).*dWf1(j,:);
            end
            dWc1  = dWf1(1,1:M) + dWf1(2,1:M); % Construct coarse brownian increments from the fine ones
            dWc2  =  dWf2(1,1:M) + dWf2(2,1:M);
            Dc = 1 + r*hc + sqrt(max(vc,0)).*dWc2;
            dRhoc = dRhoc.*Dc + Xc.*hc;
            dX0c = dX0c.*Dc;
            Xc = Xc.*Dc; 
            vc = vc + kappa*(theta-max(vc,0))*hc + xi*sqrt(max(vc,0)).*dWc1;
        end
    end
elseif model == 3  %4D GBM, euler-maruyama scheme
    d = 4;
    alf = ones(d,1)/d;
    Xf = K*ones(d,M);
    Xc = Xf;
    rho = 0.5;
    cor = eye(d) + rho*(ones(d)-eye(d));
    L = chol(cor,'lower');
    if l==0
        dWf = sqrt(hf)*L*randn(d,M);
        Xf  = Xf + r*Xf*hf + sig*Xf.*dWf;
    else
        for n = 1:nc
            for m = 1:2
                dWf = sqrt(hf)*L*randn(d,M);
                Xf  = Xf + r*Xf*hf + sig*Xf.*dWf;
                if m==1
                    dWc = dWf;
                else
                    dWc = dWc + dWf;
                end
            end
            Xc  = Xc + r*Xc*hc + sig*Xc.*dWc;
        end
    end  
end

%Compute payoffs

if option == 1  % Eurooean call
    
    ff = exp(-r*T)*max(alf'*Xf-K,0);
    fc = exp(-r*T)*max(alf'*Xc-K,0);
    
elseif option == 2 % Digital call
    
    ff = exp(-r*T)*1/2*(sign(alf'*Xf-K)+1);
    fc = exp(-r*T)*1/2*(sign(alf'*Xc-K)+1);
    
elseif option == 4 && model ~=3 % Delta european call
    
    ff = exp(-r*T)*1/2*(sign(Xf-K)+1).*dX0f;
    fc = exp(-r*T)*1/2*(sign(Xc-K)+1).*dX0c;
    
elseif option == 5 && model == 1  % GBM Vega european call
    
    ff = exp(-r*T)*1/2*(sign(Xf-K)+1).*dSigf;
    fc = exp(-r*T)*1/2*(sign(Xc-K)+1).*dSigc;
    
elseif option == 5 && model == 2 % heston Rho european call
    
    ff = exp(-r*T).*(1/2*(sign(Xf-K)+1).*dRhof-T*max(Xf-K,0));
    fc = exp(-r*T).*(1/2*(sign(Xc-K)+1).*dRhoc-T*max(Xc-K,0));
    
end

if l==0
    dF = ff;
else
    dF = ff-fc;
end

% return sum of differences, squared diff ect.

sums(1) = sum(dF);
sums(2) = sum(dF.^2);
sums(3) = sum(dF.^3);
sums(4) = sum(dF.^4);
sums(5) = sum(ff);
sums(6) = sum(ff.^2);

cost = M*nf;

end


