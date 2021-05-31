function [x,w] = GaussLaguerre(n, alpha, a, b)
%   J. Burkardt, http://people.sc.fsu.edu/~jburkardt/m_src/gen_laguerre_rule/gen_laguerre_rule.html 
%
%   Computes abscissas (x) and weights (w) of Generalized Gauss-Laguerre quadrature
%   on interval (a,+oo) for integral (x-a)^alpha * exp(-b(x-a) * f(x) dx.
%           
% Golub-Welsch technique
aj = 2*(1:n)-1 + alpha;
k  = 1:n-1;
bj = sqrt(k.*(k + alpha));
J  = diag(aj,0) + diag(bj,1) + diag(bj,-1);

[V,D] = eig(J);
[x,k] = sort(diag(D));

mu = gamma(alpha + 1.0);
w  = mu * subsref(V,substruct('()',{1,k}))'.^2; % V(1,k)

% update abcissas and weights for non-trivial arguments.
x = a + 1/b*x;
w = w.*(1/b)^(alpha+1);

end