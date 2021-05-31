%Heuristic estimate of Nq 
function [Nq] = suboptimal_Nq(eps,l,model,option,TOL_Newton)

if model == 1 || model == 3
    M = 50000;
else
    M = 100000;
end

T=1;
nf=2^l;
hf=T/nf;

If=zeros(1,M);
Ic=zeros(1,M);

res = 1;

if model == 1
    Nq = 2;
    while res > 1/(sqrt(3))*eps
        parfor i=1:M
            r = randn(nf,1);
            dWf = bb(r,T)';
            If(i) = estimator_GBM(dWf,hf,Nq,TOL_Newton,option);
            Ic(i) = estimator_GBM(dWf,hf,Nq+1,TOL_Newton,option);
        end
        Nq = Nq + 1;
        res = max(abs(If-Ic));
    end
elseif model == 2
    Nq = 10;
    while res > 1/(sqrt(3))*eps
        parfor i=1:M   
            r = randn(nf,2);
            dWf = bb(r,T)';
            If(i) = estimator_heston(dWf,hf,Nq,TOL_Newton,option);
            Ic(i) = estimator_heston(dWf,hf,Nq+1,TOL_Newton,option);
        end
        Nq = Nq + 1;
        res = max(abs(If-Ic));
    end
elseif model == 3
    Nq = 10;
    d = 4;
    rho = 0.5;
    cor = eye(4) + rho*(ones(4)-eye(4));
    L = chol(cor,'lower');
    A = getRotation(d);
    Ainv = inv(A);
    while res > 1/(sqrt(3))*eps
        parfor i=1:M
            r = randn(nf,d);
            dWf = bb(r,T)';
            If(i) = estimator_4D_GBM(dWf,hf,A,Ainv,L,Nq,TOL_Newton,option);
            Ic(i) = estimator_4D_GBM(dWf,hf,A,Ainv,L,Nq+1,TOL_Newton,option);
        end
        Nq = Nq + 1;
        res = max(abs(If-Ic));
    end
    
end

end
