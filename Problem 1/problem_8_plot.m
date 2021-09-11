MAXREL = []; %here i will store the maximum values of log10(relative error)
for ii = 1:7  %loop through files
    
    %importing data
    fname = ['solutions_7_10^', int2str(ii), '.csv' ];
    data = readmatrix(fname) ; 
    
    x = data(:,1) ;
    v = data(:,2) ;
    u = data(:,3) ;
    
    if ii == 1 %creating figure for abs err
        figabs =  figure('color' , [1 1 1]);
        axs_abs = axes();    
        xlabel('x_i')
        ylabel('log_{10}(\Delta_i)' ) 
        hold on
        grid on
        title('Absolute error as function of x_i');
    end
    
    if ii<=3 %plotting abserr as funct. of x, only for n up to 10^3
    plot(axs_abs, x, log10( abs(u-v) ) ) ; 
    end
    
    if ii == 1 %same as above, for rel err
        figrel = figure('color' , [1 1 1]);
        axs_rel = axes();
        xlabel('x_i')
        ylabel('log_{10}(\epsilon_i)' ) 
        title('Relative error as function of x_i');
        hold on
        grid on
    end
    A = log10(abs((u-v)./u)); %saving the maximum
    MAXREL = [ MAXREL, max(A) ];
    if ii<=3
    plot(axs_rel, x, A ) ; 
    end
end
figure('color' , [ 1 1 1] );
plot(1:7, MAXREL, 'marker', '.', 'MarkerSize', 8, 'linestyle', 'none') ; %plotting maximum as funct. of log10(n)
grid on;



