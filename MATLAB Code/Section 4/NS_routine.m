function NS_routine(model,option)
%Routine for computing sample smoothed call option prices, densities and sensitivies for
%GBM and Heston model dynamics. The values for Nq have been predetermined using suboptimal_Nq.m for these parameters.
%
% Model = 1 (GBM model)
% Model = 2 (Heston model)
% Model = 3 (4d Basket GBM model, only option 1 and 2)
%
% Option 1 = European Call price
% Option 2 = Digital  Call price
% Option 3 = Density of model dynamics
% Option 4 = Delta (European Call)
% Option 5 = Vega/Rho (European Call)
% Option 6 = Delta (Digital Call)
% Option 7 = Vega/Rho (Digital Call)
%

T = 1; r = 0.05; sig = 0.2; rho = -0.5; v0 = 0.04; theta = 0.0025; xi = 0.1; kappa = 1.0;

%%%%%%%%%%%%%%% GBM Model %%%%%%%%%%%%%%%%%
if model == 1 && option == 1 %GBM European Call
    X0 = 100; K = 100; eps = 0.01; Nq = 10; TOL_Newton = 1e-3;
    ref = european_call(r,sig,T,X0,K,'value');
    out = 'GBM European Call';
    
elseif model == 1 && option == 2 %GBM Digital Call
    X0 = 100; K = 100; eps = 0.001; Nq = 8; TOL_Newton = 1e-3;
    ref = digital_call(r,sig,T,X0,K,'value');
    out = 'GBM Digital Call';
    
elseif model == 1 && option == 3 %GBM Density
    X0 = 1; K = 1; eps = 0.001; Nq = 1;
    den = @(x) 1./(sqrt(2*pi*T).*sig*x).*exp(-1/2*(log(x)-log(X0)-(r-sig^2/2)*T).^2./(sig^2*T));
    ref = den(K); TOL_Newton = 1e-3;
    out = 'GBM Density';
    
elseif model == 1 && option == 4 %GBM Delta European Call
    X0 = 100; K = 100; eps = 0.001; Nq = 10; TOL_Newton = 1e-3;
    ref = european_call(r,sig,T,X0,K,'delta');
    out = 'GBM Delta European Call';
    
elseif model == 1 && option == 5 %GBM Vega European call
    X0 = 100; K = 100; eps = 0.1; Nq = 12; TOL_Newton = 1e-3;
    ref = european_call(r,sig,T,X0,K,'vega');
    out = 'GBM Vega European call';
    
elseif model == 1 && option == 6 %GBM Delta digital call
    X0 = 100; K = 100; eps = 0.0001; Nq = 1; TOL_Newton =1e-3;
    ref = digital_call(r,sig,T,X0,K,'delta');
    out = 'GBM Delta digital call';
    
elseif model == 1 && option == 7 %GBM Vega digital call
    X0 = 100; K = 100; eps = 0.01; Nq = 1; TOL_Newton = 1e-3;
    ref = digital_call(r,sig,T,X0,K,'vega');
    out = 'GBM Vega digital call';
    
    
    %%%%%%%%%%%% Heston Model %%%%%%%%%%%%
elseif model == 2 && option == 1 %Heston European Call
    X0 = 100; K = 100; eps = 0.05; Nq = 28; TOL_Newton = 1e-3;
    ref = heston(r,kappa,theta,xi,rho,X0,K,v0,T,'Euro call');
    out = 'Heston European Call';
    
elseif model == 2 && option == 2 %Heston Digital Call
    X0 = 100; K = 100; eps = 0.001; Nq = 22; TOL_Newton = 1e-3;
    ref = heston(r,kappa,theta,xi,rho,X0,K,v0,T, 'digi call');
    out = 'Heston Digital Call';
    
elseif model == 2 && option == 3 % Heston Density
    X0 = 1; K = 1; eps = 0.01; Nq = 1; TOL_Newton = 1e-3;
    ref = heston(r,kappa,theta,xi,rho,X0,K,v0,T, 'density');
    out = 'Heston density';
    
elseif model == 2 && option == 4 % Heston Delta European call
    X0 = 100; K = 100; eps = 0.01; Nq = 22; TOL_Newton = 1e-3;
    ref = heston(r,kappa,theta,xi,rho,X0,K,v0,T,'Euro delta');
    out = 'Heston Delta European call';
    
elseif model == 2 && option == 5 % Heston Rho European call
    X0 = 100; K = 100; eps = 0.1; Nq = 28; TOL_Newton = 1e-3;
    ref = heston(r,kappa,theta,xi,rho,X0,K,v0,T,'Euro rho');
    out = 'Heston Rho European call';
    
elseif model == 2 && option == 6 % Heston Delta Digital call
    X0 = 100; K = 100; eps = 0.001; Nq = 1; TOL_Newton = 1e-3;
    ref = heston(r,kappa,theta,xi,rho,X0,K,v0,T,'digi delta');
    out = 'Heston Delta Digital call';
    
elseif model == 2 && option == 7 % Heston Rho Digital call
    X0 = 100; K = 100; eps = 0.01; Nq = 1; TOL_Newton = 1e-3;
    ref = heston(r,kappa,theta,xi,rho,X0,K,v0,T,'digi rho');
    out = 'Heston Rho digital call';
    
    
    %%%%%%%%%%%%%%%% 4D GBM Model %%%%%%%%%%%%%%%%%
elseif model == 3 && option == 1 %GBM European Call
    d = 4; rho = 0.5; X0 = 100; K = 100; eps = 0.05; Nq = 15; TOL_Newton = 1e-3;
    ref = basket_call(r,sig,T,X0,K,rho,d,'Euro call');
    out = '4D GBM European Call';
    
elseif model == 3 && option == 2 %GBM Digital Call
    d = 4; rho = 0.5; X0 = 100; K = 100; eps = 0.001; Nq = 13; TOL_Newton = 1e-3;
    ref = basket_call(r,sig,T,X0,K,rho,d,'digi call');
    out = '4D GBM Digital Call';
    
end

M0    = 1000;   % initial samples on coarse levels
Lmin  = 2;      % minimum refinement level
Lmax  = 12;     % maximum refinement level

[res,~,~ ] = mlmc_smooth(@estimator_l,M0,Nq,TOL_Newton,model,option,eps,Lmin,Lmax,0,0,0); %Smoothed MLMC
%[res,~,~] = mlmc(@level_estimator,M0,model,option,eps,Lmin,Lmax,0,0,0);                  %standard MLMC
%
%RMSE check
RMSE = sqrt((res-ref)^2);
fprintf('\n')
disp(['Model:               ',out])
disp(['MLMC smoothed value: ',num2str(res)])
disp(['Reference value:     ',num2str(ref)])
disp(['TOL_Newton:          ',num2str(TOL_Newton)])
disp(['Nq:                  ',num2str(Nq)])
disp(['RMSE:                ',num2str(eps)])
disp(['RMSE check:          ',num2str(RMSE)])
fprintf('\n')
end


