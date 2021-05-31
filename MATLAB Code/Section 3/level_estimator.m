function [ sums, cost] = level_estimator(l,M,option)

T = 1; X0 = 100; K = 100; r = 0.05; sig = 0.2;

sums(1:6) = 0;

nf = 2^l;
nc = nf/2;
dtf = T/nf;
dtc = T/nc;

Xf = X0*ones(1,M);
Xc = Xf;

dX0f = ones(1,M);
dX0c = dX0f;

dSigf = zeros(1,M);
dSigc = dSigf;

if option == 3 || option == 4 || option == 5 || option == 6 || option == 8
    nc = nc - 1;
end

%Gaussian CDF and PDF
Phi = @(x)(1+erf(x/sqrt(2)))/2;
phi = @(x) 1/sqrt(2*pi).*exp(-1/2*x.^2);

if l==0 && option ~= 3 && option ~= 4 && option ~= 5 && option ~= 6 && option ~= 8 % Just one timestep  
    dWf = sqrt(dtf)*randn(1,M);
    Mil_f = 1 + r*dtf + sig.*dWf + 1/2*sig^2*(dWf.^2-dtf);
    dSigf = dSigf.*Mil_f + Xf.*(dWf + sig*(dWf.^2-dtf));
    dX0f = dX0f.*Mil_f; 
    Xf = Xf.*Mil_f;  
else
    for k = 1:nc
        dWf = sqrt(dtf)*randn(2,M);
        for j = 1:2
            Mil_f = 1 + r*dtf + sig.*dWf(j,:) + 1/2*sig^2*(dWf(j,:).^2-dtf);
            dSigf = dSigf.*Mil_f + Xf.*(dWf(j,:) + sig*(dWf(j,:).^2-dtf));
            dX0f = dX0f.*Mil_f;
            Xf = Xf.*Mil_f;  
        end
        dWc = dWf(1,:) + dWf(2,:);  % Construct coarse brownian increments from the fine ones
        Mil_c = 1 + r*dtc + sig.*dWc + 1/2*sig^2*(dWc.^2-dtc);
        dSigc = dSigc.*Mil_c + Xc.*(dWc + sig*(dWc.^2-dtc));
        dX0c = dX0c.*Mil_c;
        Xc = Xc.*Mil_c;
    end
end

%Compute payoffs

if option == 1 % Delta european call
    
    ff = 1/2*(sign(Xf-K)+1).*dX0f;
    fc = 1/2*(sign(Xc-K)+1).*dX0c;
    
elseif option == 2 % Vega european call
    
    ff = 1/2*(sign(Xf-K)+1).*dSigf;
    fc = 1/2*(sign(Xc-K)+1).*dSigc;
    
elseif option == 3 || option == 5  % Smoothed Delta european and digital call
    
    if(l==0)
        af = (1+r*dtf).*Xf;
        bf = sig.*Xf*sqrt(dtf);

        daf = (1+r*dtf).*dX0f;
        dbf = sig*dX0f.*sqrt(dtf);
       
        if option == 3
            ff = daf.*Phi((af-K)./bf) + dbf.*phi((af-K)./bf);
        elseif option == 5
            ff = ((daf.*bf-(af-K).*dbf)./(bf.^2)).*phi((af-K)./bf);
        end
        fc  = ff;
    else
        dWf = sqrt(dtf)*randn(1,M);
        Mil = 1 + r*dtf + sig.*dWf + 1/2*sig^2*(dWf.^2-dtf);
        dX0f = dX0f.*Mil;
        Xf = Xf.*Mil;
        
        af = (1+r*dtf).*Xf;
        bf = sig.*Xf*sqrt(dtf);
        daf = (1+r*dtf).*dX0f;
        dbf = sig*dX0f.*sqrt(dtf);
             
        ac = (1+r*dtc + sig.*dWf).*Xc;
        bc = sig.*Xc*sqrt(dtf);
        dac = (1+r*dtc + sig.*dWf).*dX0c;
        dbc = sig*dX0c.*sqrt(dtf);
        
        if option == 3
            ff = daf.*Phi((af-K)./bf) + dbf.*phi((af-K)./bf);
            fc = dac.*Phi((ac-K)./bc) + dbc.*phi((ac-K)./bc); 
        elseif option == 5
            ff = ((daf.*bf-(af-K).*dbf)./(bf.^2)).*phi((af-K)./bf);
            fc = ((dac.*bc-(ac-K).*dbc)./(bc.^2)).*phi((ac-K)./bc);
        end
     
    end
    
elseif option == 4 || option == 6  % Smoothed Vega european and digital call
    
    if(l==0)
        af = (1+r*dtf).*Xf;
        bf = sig.*Xf*sqrt(dtf);
        daf = (1+r*dtf).*dSigf ;
        dbf = sig*dSigf.*sqrt(dtf) + Xf.*sqrt(dtf);
        
        if option == 4
            ff = daf.*Phi((af-K)./bf) + dbf.*phi((af-K)./bf);
        elseif option == 6
            ff = ((daf.*bf-(af-K).*dbf)./(bf.^2)).*phi((af-K)./bf);
        end
        fc  = ff;
    else
        dWf = sqrt(dtf)*randn(1,M);
        Mil = 1 + r*dtf + sig.*dWf + 1/2*sig^2*(dWf.^2-dtf);
        
        dSigf = dSigf.*Mil + Xf.*(dWf + sig*(dWf.^2-dtf));
        Xf = Xf.*Mil;
        
        af = (1+r*dtf).*Xf;
        bf = sig.*Xf*sqrt(dtf);
        daf = (1+r*dtf).*dSigf;
        dbf = sig*dSigf.*sqrt(dtf) + Xf.*sqrt(dtf);
        
        ac = (1+r*dtc + sig.*dWf).*Xc;
        bc = sig.*Xc*sqrt(dtc/2);
        dac = (1+r*dtc + sig*dWf).*dSigc + Xc.*dWf;
        dbc = sig*dSigc.*sqrt(dtc/2) + Xc*sqrt(dtc/2);
        
        if option == 4
            ff = daf.*Phi((af-K)./bf) + dbf.*phi((af-K)./bf);
            fc = dac.*Phi((ac-K)./bc) + dbc.*phi((ac-K)./bc);
        elseif option == 6
            ff = ((daf.*bf-(af-K).*dbf)./(bf.^2)).*phi((af-K)./bf);
            fc = ((dac.*bc-(ac-K).*dbc)./(bc.^2)).*phi((ac-K)./bc);
        end
    end
    
elseif option == 7 %digital call
    
    ff = 1/2*(sign(Xf-K)+1);
    fc = 1/2*(sign(Xc-K)+1);
    
elseif option == 8 %smooted digital call
    if(l==0)
        ff  = Phi((Xf+r*Xf*dtf-K)./(sig*Xf*sqrt(dtf)));
        fc  = ff;
    else
        dWf = sqrt(dtf)*randn(1,M);
        Xf  = Xf + r*Xf*dtf + sig*Xf.*dWf + 1/2*sig^2*Xf.*(dWf.^2-dtf);
        ff  = Phi((Xf+r*Xf*dtf-K)./(sig*Xf*sqrt(dtf)));
        fc  = Phi((Xc+r*Xc*dtc+sig*Xc.*dWf-K)./(sig*Xc*sqrt(dtc/2)));
    end
   
end

ff = exp(-r*T)*ff;
fc = exp(-r*T)*fc;

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


