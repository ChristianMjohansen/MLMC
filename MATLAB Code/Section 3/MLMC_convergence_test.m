function MLMC_convergence_test

close all;
addpath('..');
Lmin  = 2;      % minimum refinement level
Lmax  = 10;     % maximum refinement level

T = 1; X0 = 100; K = 100; r = 0.05; sig = 0.2;

M    = 10^7;    % samples for convergence tests
L    = 8;       % levels for convergence tests

for opt = 1 % Choose option
    if opt == 1
        fprintf(1,'\n ---- European call Delta ---- \n');
        M0   = 5000;                                % Initial samples 
        Eps  = [ 1e-4 3e-4 5e-4 7e-4 ];             % RMSE tolerance
        ref = european_call(r,sig,T,X0,K,'delta');  % Reference value
        nvert = 1;
        filename = ['Greek_' num2str(opt)];
    elseif opt == 2
        fprintf(1,'\n ---- European call Vega ---- \n');
        M0   = 5000;        
        Eps  = [ 2e-2 3e-2 4e-2 5e-2 ];
        ref  = european_call(r,sig,T,X0,K,'vega');
        nvert = 1;
        filename = ['Greek_' num2str(opt)];
    elseif opt == 3
        fprintf(1,'\n ---- European call Delta (Smoothed)---- \n');
        M0   = 1000; 
        Eps  = [ 1e-4 3e-4 5e-4 7e-4 ];
        ref  = european_call(r,sig,T,X0,K,'delta');
        nvert = 1;
        filename = ['Greek_' num2str(opt)];
    elseif opt == 4
        fprintf(1,'\n ---- Euroean call Vega (Smoothed) ---- \n');
        M0   = 1000; 
        Eps  = [ 2e-2 3e-2 4e-2 5e-2 ];
        ref  = european_call(r,sig,T,X0,K,'vega');
        nvert = 1;
        filename = ['Greek_' num2str(opt)];
    elseif opt == 5
        fprintf(1,'\n ---- Digital call Delta (Smoothed) ---- \n');
        M0   = 1000; 
        Eps  = [ 1e-5 3e-5 5e-5 7e-5 ];
        ref  = digital_call(r,sig,T,X0,K,'delta');
        nvert = 1;
        filename = ['Greek_' num2str(opt)];
    elseif opt == 6
        fprintf(1,'\n ---- Digital call Vega (Smoothed) ---- \n');
        M0   = 1000; 
        Eps  = [ 1e-4 3e-4 5e-4 7e-4 ];
        ref  = digital_call(r,sig,T,X0,K,'vega'); 
        nvert = 1;
        filename = ['Greek_' num2str(opt)];
    elseif opt == 7
        fprintf(1,'\n ---- Digital call option price  ---- \n');
        M0    = 5000;
        Eps  = [ 5e-5 1e-4 2e-4 4e-4 ];
        ref = digital_call(r,sig,T,X0,K,'value');
        nvert = 2;
        filename = 'Digital_call';
    elseif opt == 8
        fprintf(1,'\n ---- Digital call option price (Smoothed) ---- \n');
        M0    = 1000;
        Eps  = [ 5e-5 1e-4 2e-4 4e-4 ];
        ref = digital_call(r,sig,T,X0,K,'value');
        nvert = 2;
        filename = 'Smoothed_digital_call';
    end
    
    fp = fopen(['results/' filename '.txt'],'w');
    mlmc_test(@level_estimator, M, L, M0,opt, Eps, Lmin, Lmax, fp);
    fclose(fp);
    
    fprintf(1,'\n Reference value: %f \n\n',ref);
    
    %
    % plot results
    %
    mlmc_plot(filename, nvert);
    print(gcf,'-depsc','-painters',['results/' filename '.eps'])
end

end