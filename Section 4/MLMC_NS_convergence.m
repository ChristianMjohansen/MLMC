% MLMC convergence tests for standard MLMC and MLMC with numerical
% smoothing
function MLMC_NS_convergence

close all;

addpath('..');

Lmin  = 2;      % minimum refinement level
Lmax  = 12;     % maximum refinement level
TOL_Newton = 1e-3;

for opt = 13 % select option
    if opt == 1
        fprintf(1,'\n ---- GBM digital call ---- \n');
        filename = 'GBM_digital_call';
        model = 1; option = 2;
        M     = 10^6;                   % samples for convergence tests
        M2    = 10^7;                   % samples for convergence tests without smoothing
        M0    = 1000;                   % initial samples
        L     = 8;                      % levels for convergence tests
        %Eps   = [ 5e-4 1e-3 5e-3 1e-2];
        %Nq    = [ 14 11 8 6 ]; % pre-calculated with suboptimal_Nq.m
        Eps   =  5e-4;
        Nq    =  14;
    elseif opt == 2
        fprintf(1,'\n ---- Heston digital call ---- \n');
        filename = 'Heston_digital_call';
        model = 2; option = 2;
        M     = 10^7;
        M2    = 10^8;
        M0    = 1000;
        L     = 8;
        % Eps   = [ 5e-4 1e-3 5e-3 1e-2 ];
        % Nq    = [ 32 28 21 19 ];
        Eps   =  5e-4;
        Nq    =  32;
    elseif opt == 3
        fprintf(1,'\n ---- 4D GBM digital call ---- \n');
        filename = '4D_GBM_digital_call';
        model = 3; option = 2;
        M     = 10^6;
        M2    = 10^8;
        M0    = 1000;
        L     = 8;
        %Eps   = [ 5e-4 1e-3 5e-3 1e-2 ];
        %Nq    = [ 20 17 13 11 ];
        Eps   =  5e-4;
        Nq    =  20;
    elseif opt == 4
        fprintf(1,'\n ---- GBM Density ---- \n');
        filename = 'GBM_density';
        model = 1; option = 3;
        M     = 10^6;
        M0    = 1000;
        L     = 9;
        Eps   =  1e-3;
        Nq    =  1;
    elseif opt == 5
        fprintf(1,'\n ---- Heston Density  ---- \n');
        filename = 'Heston_density';
        model = 2; option = 3;
        M     = 10^6;
        M0    = 1000;
        L     = 8;
        %Eps   = [ 1e-3 5e-3 1e-2 ];
        Eps   = 1e-3;
        Nq    =  1;
        
%%%%%%%%%%%% GBM Greeks %%%%%%%%%%%%
    elseif opt == 6
        fprintf(1,'\n ---- GBM European call Delta  ---- \n');
        filename = 'GBM_euro_delta';
        model = 1; option = 4;
        M     = 10^6;
        M2    = 10^7;
        M0    = 1000;
        L     = 8;
        Eps   = 1e-3;
        Nq    = 10;
    elseif opt == 7
        fprintf(1,'\n ---- GBM European call Vega  ---- \n');
        filename = 'GBM_euro_vega';
        model = 1; option = 5;
        M     = 10^6;
        M0    = 1000;
        M2    = 10^7;
        L     = 8;
        Eps   = 0.1;
        Nq    = 12;
    elseif opt == 8
        fprintf(1,'\n ---- GBM Digital call Delta   ---- \n');
        filename = 'GBM_dig_delta';
        model = 1; option = 6;
        M     = 10^6;
        M0    = 1000;
        L     = 9;
        Eps   = 1e-3;
        Nq    = 1;
    elseif opt == 9
        fprintf(1,'\n ---- GBM Digital call Vega  ---- \n');
        filename = 'GBM_dig_vega';
        model = 1; option = 7;
        M     = 10^6;
        M0    = 1000;
        L     = 8;
        Eps   = 1e-2;
        Nq    = 1;
        
%%%%%%%%%%%% Heston Greeks  %%%%%%%%%%%%
    elseif opt == 10
        fprintf(1,'\n ---- Heston European call Delta  ---- \n');
        filename = 'Heston_euro_delta';
        model = 2; option = 4;
        M     = 10^6;
        M2    = 10^8;
        M0    = 1000;
        L     = 8;
        Eps   = 1e-2;
        Nq    = 22;
    elseif opt == 11
        fprintf(1,'\n ---- Heston European call Rho  ---- \n');
        filename = 'Heston_euro_rho';
        model = 2; option = 5;
        M     = 10^6;
        M0    = 1000;
        M2    = 10^7;
        L     = 8;
        Eps   = 0.1;
        Nq    = 28;
    elseif opt == 12
        fprintf(1,'\n ---- Heston Digital call Delta   ---- \n');
        filename = 'Heston_dig_delta';
        model = 2; option = 6;
        M     = 10^6;
        M0    = 1000;
        L     = 9;
        Eps   = 1e-3;
        Nq    = 1;
    elseif opt == 13
        fprintf(1,'\n ---- Heston Digital call Rho  ---- \n');
        filename = 'Heston_dig_rho';
        model = 2; option = 7;
        M     = 10^6;
        M0    = 1000;
        L     = 8;
        Eps   = 1e-2;
        Nq    = 1;
    end
% Uncomment either "Smoothed MLMC" or "standard MLMC"   
    
    %Smoothed MLMC
    % fp = fopen(['Results/' filename '.txt'],'w');
    % mlmc_test(@estimator_l, M, L, M0, Nq, TOL_Newton, model,option, Eps, Lmin, Lmax, fp);
    % fclose(fp);
    
    %Standard MLMC
     fp2 = fopen(['Results/without_NS/' filename '.txt'],'w');
     mlmc_test2(@level_estimator, M2, L, M0, model, option, Eps, Lmin, Lmax, fp2 )
     fclose(fp2);
    
    %
    % plot results
    %
    
    nvert = 2; arg = 2;
    mlmc_plot(filename, nvert, arg);
    %print(gcf,'-depsc','-painters',['Results/' filename '.eps'])
    print(gcf,'-depsc','-painters',['Results/without_NS/' filename '.eps'])
    
end

end