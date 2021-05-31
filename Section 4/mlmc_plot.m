%
% utility to generate MLMC plots based on input text file mlmc_plot(filename,nvert)
%

function mlmc_plot(filename,nvert,arg,varargin)

close all;

if nargin==3
    error_bars = 0;
elseif nargin==4
    if varargin{1}=='error_bars'
        error_bars = 1;
    else
        error('invalid mlmc_plot option')
    end
else
    nargin
    error('invalid number of mlmc_plot arguments')
end

%
% read in data
%
if arg == 1
    fid = fopen(['Results/' filename '.txt'],'r');
elseif arg == 2
    fid = fopen(['Results/without_NS/' filename '.txt'],'r');
end

line = '    ';
while (length(line)<20) | (strcmp(line(1:4),'*** ')==0)
    line = [ fgetl(fid) '    ' ];
end
file_version = sscanf(line(23:30),'%f');
if isempty(file_version)
    file_version = 0.8;
end

if (file_version<0.9)
    M = 0;
    if (error_bars)
        error('cannot plot error bars -- no value of M in file');
    end
else
    while (length(line)<20) | (strcmp(line(1:9),'*** using')==0)
        line = [ fgetl(fid) '    ' ];
    end
    M = sscanf(line(14:20),'%d');
end

line = '    ';
while (length(line)<20) | (strcmp(line(1),'-')==0)
    line = [ fgetl(fid) '    ' ];
end

line = fgetl(fid);
l    = 1;
while (length(line)>10)
    data = sscanf(line,'%f');
    del1(l) = data(2);
    del2(l) = data(3);
    var1(l) = data(4);
    var2(l) = data(5);
    kur1(l) = data(6);
    chk1(l) = data(7);
    cost(l) = data(8);
    
    line = fgetl(fid);
    l    = l+1;
end

vvr1 = var1.^2 .* (kur1-1);

L = l-2;

line = '    ';
while (length(line)<20) | (strcmp(line(1),'-')==0)
    line = [ fgetl(fid) '    ' ];
end

line = fgetl(fid);
l    = 1;

while (length(line)>10)
    data = sscanf(line,'%f');
    Eps(l)       = data(1);
    mlmc_cost(l) = data(3);
    std_mlmc_cost(l)  = data(4);
    len          = length(data)-5;
    ls(1:len,l)  = 0:len-1;
    Mls(1:len,l) = data(6:end);
    
    line = fgetl(fid);
    l    = l+1;
end

%
% plot figures
%
L1 = 1:L;

figs(1) = figure;

pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.75*nvert]; set(gcf,'pos',pos);

set(0,'DefaultAxesColorOrder',[0 0 1; 0 0 1; 1 0 0]);
set(0,'DefaultAxesLineStyleOrder','-*|:*|--o|--x|--d|--s')
%set(0,'DefaultAxesLineStyleOrder','-*|--*')

if nvert == 2
    subplot(nvert,2,1)
    plot(0:L,log2(var2),'-+',1:L,log2(var1(2:end)),'--*',L1,log2(2.^(-L1))-4.5,'--o')
    xlabel('level $\ell$','Interpreter','latex');
    ylabel('$\log_2$ variance','Interpreter','latex');
    
    current = axis; axis([ 0 L current(3:4) ]);
    legend({'$P_\ell$','$P_\ell\!-\! P_{\ell-1}$','$2^{-\ell}$'}, ...
        'Interpreter','latex','Location','SouthWest')
    grid on
    hold on;
    
    if error_bars
        plot([1:L; 1:L],[log2(max(abs(var1(2:end))-3*sqrt(vvr1(2:end)/M),1e-10)); ...
            log2(    abs(var1(2:end))+3*sqrt(vvr1(2:end)/M))],'-r.')
    end
    
    subplot(nvert,2,2)
    plot(0:L,log2(abs(del2)),'-+',1:L,log2(abs(del1(2:end))),'--*',L1,log2(2.^(-L1))-3.2,'--o')
    xlabel('level $\ell$','Interpreter','latex');
    ylabel('$\log_2 |\mbox{mean}|$','Interpreter','latex');
    current = axis; axis([ 0 L current(3:4) ]);
    legend({'$P_\ell$','$P_\ell\!-\! P_{\ell-1}$','$2^{-\ell}$'}, ...
        'Interpreter','latex','Location','Northeast')
    grid on
    hold on;
    
    %plot([0:L; 0:L],[log2(max(abs(del2)-3*sqrt(var2/M),1e-10)); ...
    %                 log2(    abs(del2)+3*sqrt(var2/M))],'-r.')
    if error_bars
        plot([1:L; 1:L],[log2(max(abs(del1(2:end))-3*sqrt(var1(2:end)/M),1e-10)); ...
            log2(    abs(del1(2:end))+3*sqrt(var1(2:end)/M))],'-r.')
    end
    set(0,'DefaultAxesColorOrder',[0 0 1; 1 0 0]);
    
    current = axis; axis([ 0 L  current(3:4) ]);
    
    subplot(nvert,2,3)
    plot(0:L,log2(cost(1:end)),'--*')
    grid on
    hold on
    plot(L1,log2(2.^(L1))+1.4,'--o')
    xlabel('level $\ell$','Interpreter','latex');
    ylabel('$\log_2$ cost per level','Interpreter','latex');
    legend({'$\log_2$ cost per level','$2^{\ell}$'}, ...
        'Interpreter','latex','Location','NorthWest')
    current = axis; axis([ 0 L  current(3:4) ]);
    
    subplot(nvert,2,4)
    semilogy(1:L,kur1(2:end),'--*')
    xlabel('level $\ell$','Interpreter','latex');
    ylabel('kurtosis');
    current = axis; axis([ 0 L current(3:4) ]);
    grid on
    
%     set(0,'DefaultAxesColorOrder',[0 0 1]);
%     subplot(nvert,2,5)
%     semilogy(ls, Mls)
%     xlabel('level $\ell$','Interpreter','latex');
%     ylabel('$N_\ell$','Interpreter','latex');
%     current = axis; axis([ 0 size(Mls,1)-1 current(3:4) ]);
%     for i=1:length(Eps)
%         labels{i} = num2str(Eps(i));
%     end
%     legend(labels,'Location','NorthEast')
%     
%     grid on
%     
%     set(0,'DefaultAxesColorOrder',[0 0 1]);
%     subplot(nvert,2,6)
%     plot(1:L,chk1(2:end),'--*')
%     xlabel('level l'); ylabel('consistency check');
%     set(0,'DefaultAxesLineStyleOrder','-*|:*')
%     %set(0,'DefaultAxesLineStyleOrder','-*|--*')
%     
end

if nvert==1
    pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.75]; set(gcf,'pos',pos);
    
    loglog(Eps, mlmc_cost(:)','--ro',Eps,std_mlmc_cost(:)','--bo' , Eps, Eps.^(-2).*log(Eps).^2*1/2,'-k',  Eps, Eps.^(-2.5)*2,'-g' )
    xlabel('accuracy $\epsilon$','Interpreter','latex');
    ylabel('Expected Cost','Interpreter','latex');
    current = axis; axis([ Eps(1) Eps(end) current(3:4) ]);
    legend('MLMC + smoothing', 'MLMC without smoothing', '$\epsilon^{-2}\log(\epsilon)^2$', '$\epsilon^{-2.5}$', 'Interpreter','latex')
    grid on
    
    set(0,'DefaultAxesLineStyleOrder',':o|:x|:d|:*|:s');
    %set(0,'DefaultAxesLineStyleOrder','--o|--x|--d|--*|--s');
end


