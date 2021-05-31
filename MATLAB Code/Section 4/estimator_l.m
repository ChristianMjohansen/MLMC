% Level estimator for MLMC with numerical smoothing
function [sums, cost] = estimator_l(l,M,Nq,TOL_Newton,model,option)

sums(1:6) = 0;

If = zeros(1,M);
Ic = zeros(1,M);

nf = 2^l;
nc = nf/2;

T = 1;
hf = T/nf;
hc = T/nc;

if model == 1
    d = 1;
    if l==0
        parfor i=1:M
            r = randn(nf,d);
            dWf = bb(r,T)'; % Brownian Bridge construction
            If(i) = estimator_GBM(dWf,hf,Nq,TOL_Newton,option);
        end
    else
        parfor i=1:M
            r = randn(nf,d);
            dWf =  bb(r,T)';
            dWc = dWf(:,1:2:end) + dWf(:,2:2:end);
            If(i) = estimator_GBM(dWf,hf,Nq,TOL_Newton,option);
            Ic(i) = estimator_GBM(dWc,hc,Nq,TOL_Newton,option);
        end
    end
elseif model == 2
    d = 2;
    if l==0
        parfor i=1:M
            r = randn(nf,d);
            dWf = bb(r,T)';
            If(i) = estimator_heston(dWf,hf,Nq,TOL_Newton,option);
        end
    else
        parfor i=1:M
            r = randn(nf,d);
            dWf =  bb(r,T)';
            dWc = dWf(:,1:2:end) + dWf(:,2:2:end);
            If(i) = estimator_heston(dWf,hf,Nq,TOL_Newton,option);
            Ic(i) = estimator_heston(dWc,hc,Nq,TOL_Newton,option);
        end
    end
elseif model == 3
    d = 4;
    rho = 0.5;
    cor = eye(d) + rho*(ones(d)-eye(d));
    L = chol(cor,'lower');
    A = getRotation(d);
    Ainv = inv(A);
    if l==0
        parfor i=1:M
            r = randn(nf,d);
            dWf = bb(r,T)';
            If(i) = estimator_4D_GBM(dWf,hf,A,Ainv,L,Nq,TOL_Newton,option);
        end
    else
        parfor i=1:M
            r = randn(nf,d);
            dWf =  bb(r,T)';
            dWc = dWf(:,1:2:end) + dWf(:,2:2:end);
            If(i) = estimator_4D_GBM(dWf,hf,A,Ainv,L,Nq,TOL_Newton,option);
            Ic(i) = estimator_4D_GBM(dWc,hc,A,Ainv,L,Nq,TOL_Newton,option);
        end
    end
end

if l==0
    dI = If;
else
    dI = If - Ic;
end

sums(1) = sum(dI);
sums(2) = sum(dI.^2);
sums(3) = sum(dI.^3);
sums(4) = sum(dI.^4);
sums(5) = sum(If);
sums(6) = sum(If.^2);

%cost = M*nf;   % cost
cost = M*nf*Nq*log(1/TOL_Newton);
end