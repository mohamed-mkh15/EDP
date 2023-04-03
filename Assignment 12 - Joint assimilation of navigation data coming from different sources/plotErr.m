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
    plot(6:N,ErrAvg(4:end))
    hold on
    plot((m+4):N,ErrAvgE(4:end),'r')

    plot(P,'color','m','linewidth',1.5)
    title(sprintf('%s errors',s))
    legend('filteration True Estimation error','Extrapolation True estimation error','Standard deviation')
    set(gcf,'position',[0,0,900,800]);
    
    xlabel('points')
    ylabel('value')
    grid on
end

