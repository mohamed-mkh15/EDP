function [  ] = plotErr( Err,ErrE,P, s )

%PlotErr This function is to calculate and plot average errors
%   This function is to calculate then plot the true estimation error
%   for Kalman filter filtered data and extrapolated data
%   Output arguments: NaN, only Graphs and visulaization
%   Input arguments:
%   Err     True estimation error matrix (over iterations)
%   ErrE    Extrapolation True estimation error matrix (over iterations)
%   P       Calculation error vector (standard deviation calculated from 
%           covariance matrix
%   s       String containing the variable name

    [M,N]=size(Err);     %number of iterations and points
    [~,NE]=size(ErrE);   %number of extrapolation error points
    m=N-NE;              %extrapolation steps
    
    ErrAvg=zeros(1,N-2);   %Average true estimation error
    for j=1:N-2;
        ErrAvg(j)=sqrt((1/(M-1))*sum(Err(:,j+2)));
    end
    ErrAvgE=zeros(1,NE);  %Average extrapolation true estimation error
    for j=1:NE;
        ErrAvgE(j)=sqrt((1/(M-1))*sum(ErrE(:,j)));
    end
    
    %plotting
    plot(3:N,ErrAvg)
    hold on
    plot((m+1):N,ErrAvgE,'r')
    if isnan(P)
        title(sprintf('%s True Estimation Error vs. Extrapolation error',s))
        legend('filteration True Estimation error','Extrapolation True estimation error')
    else
        plot(P,'m')
        title(sprintf('%s True Estimation Error vs. Extrapolation error vs. calculation error',s))
        legend('filteration True Estimation error','Extrapolation True estimation error','calculation error')
    end
    xlabel('points')
    ylabel('value')
    grid on
end

