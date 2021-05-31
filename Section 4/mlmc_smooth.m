% function [P, Ml, cost] = mlmc_smooth(mlmc_l,M0,eps,Lmin,Lmax,
%                               alpha,beta,gamma)
%
% multi-level Monte Carlo estimation for numerical smoothing
%
% P     = value
% Ml    = number of samples at each level
% cost  = total cost
%
% M0    = initial number of samples         > 0
% eps   = desired accuracy (rms error)      > 0 
% Lmin  = minimum level of refinement       >= 2
% Lmax  = maximum level of refinement       >= Lmin
%
% alpha -> weak error is  O(2^{-alpha*l})
% beta  -> variance is    O(2^{-beta*l})
% gamma -> sample cost is O(2^{gamma*l})
%
% if alpha, beta, gamma are not positive, then they will be estimated
%
% mlmc_l = function for level l estimator 
%
% inputs:  l = level
%          M = number of samples
%         
% output: sums(1) = sum(Y)
%         sums(2) = sum(Y.^2)
%         where Y are iid samples with expected value:
%         E[P_0]           on level 0
%         E[P_l - P_{l-1}] on level l>0
%         cost = cost of M samples

function [P, Ml, Cl] = mlmc_smooth(mlmc_l, M0, Nq, TOL_Newton, model, option, eps, Lmin, Lmax, ...
                            alpha0, beta0, gamma0)

%
% check input parameters
%
  if (Lmin<2)
    error('error: needs Lmin >= 2');
  end

  if (Lmax<Lmin)
    error('error: needs Lmax >= Lmin');
  end

  if (M0<=0 || eps<=0)
    error('error: needs M0>0, eps>0 \n');
  end

%
% initialisation
%
  alpha = max(0, alpha0);
  beta  = max(0, beta0);
  gamma = max(0, gamma0);


  L = Lmin;

  Ml(1:L+1)       = 0;
  suml(1:2,1:L+1) = 0;
  costl(1:L+1)    = 0;
  dMl(1:L+1)      = M0;

  while sum(dMl) > 0

%
% update sample sums
%
    for l=0:L
      if dMl(l+1) > 0
        [sums, cost] = mlmc_l(l,dMl(l+1),Nq, TOL_Newton,model,option);
        Ml(l+1)     = Ml(l+1)     + dMl(l+1);
        suml(1,l+1) = suml(1,l+1) + sums(1);
        suml(2,l+1) = suml(2,l+1) + sums(2);
        costl(l+1)  = costl(l+1)  + cost;
      end
    end

%
% compute absolute average, variance and cost
%
    ml = abs(   suml(1,:)./Ml);
    Vl = max(0, suml(2,:)./Ml - ml.^2);
    Cl = costl./Ml;

%
% fix to cope with possible zero values for ml and Vl
% (can happen in some applications when there are few samples)
%
    for l = 3:L+1
      ml(l) = max(ml(l), 0.5*ml(l-1)/2^alpha);
      Vl(l) = max(Vl(l), 0.5*Vl(l-1)/2^beta);
    end

%
% use linear regression to estimate alpha, beta, gamma if not given
%
    A = repmat((1:L)',1,2).^repmat(1:-1:0,L,1);

    if alpha0 <= 0
      x     = A \ log2(ml(2:end))';
      alpha = max(0.5,-x(1));
    end

    if beta0 <= 0
      x     = A \ log2(Vl(2:end))';
      beta  = max(0.5,-x(1));
    end

    if gamma0 <= 0
      x     = A \ log2(Cl(2:end))';
      gamma = max(0.5,x(1));
    end
%
% set optimal number of additional samples
%
    Ms  = ceil( 3*sqrt(Vl./Cl)*sum(sqrt(Vl.*Cl)) /  eps^2 );
    dMl = max(0, Ms-Ml);
%
% if (almost) converged, estimate remaining error and decide 
% whether a new level is required
%
    if sum( dMl > 0.01*Ml ) == 0
      % this change copes with cases with erratic ml values
      %range = 0:min(2,L-1);
      %rem = max(ml(L+1-range) ./ 2.^(range*alpha)) / (2^alpha - 1);
      rem = ml(L+1) / (2^alpha - 1);

      if rem > eps/sqrt(3)
        if (L==Lmax)
          fprintf(1,'*** failed to achieve weak convergence *** \n');
        else
          L       = L+1;
          Vl(L+1) = Vl(L) / 2^beta;
          Cl(L+1) = Cl(L) * 2^gamma;
          Ml(L+1) = 0;
          suml(1:2,L+1) = 0;
          costl(L+1) = 0;

          Ms  = ceil(3*sqrt(Vl./Cl) * sum(sqrt(Vl.*Cl)) / eps^2);
          dMl = max(0, Ms-Ml);
        end
      end
    end
  end

%
% finally, evaluate multilevel estimator
%
  P = sum(suml(1,:)./Ml);
end