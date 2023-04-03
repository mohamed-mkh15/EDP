function [] = M_runs(MM,p)
    Err=zeros(500,200);
    ErrV=zeros(500,200);
    ErrE1=zeros(500,199);
    ErrE7=zeros(500,193);
    ErrVE1=zeros(500,199);
    ErrVE7=zeros(500,193);
    for M=1:MM
        normaldist=makedist('Normal',0,0.2);

        a=random(normaldist,200,1);

        X=zeros(1,200);
        X(1)=5;
        V=zeros(1,200);
        V(1)=1;
        T=1;
        N=length(X);
        for i=2:N
            X(i)=X(i-1)+V(i-1)*T+a(i-1)*T^2/2;
            V(i)=V(i-1)+a(i-1)*T;
        end

        normaldist2=makedist('Normal',0,20);
        eta=random(normaldist2,200,1);
        Z=zeros(1,200);
        N=length(X);
        Z(1)=X(1)+eta(1);
        for i=2:N
            zeta=rand;
            if zeta>p
                Z(i)=X(i)+eta(i);
            else
                Z(i)=NaN;
            end
        end

        phi=[1 T;0 1];
        G=[T^2/2; T];
        H=[1 0];
        sigma_a=0.2^2;
        Q=G*G'*sigma_a;

        xi=zeros(1,200);
        
        xiE1=zeros(1,200);
        xiE7=zeros(1,200);
        
        Xi=[2;0];
        P=[10000 0; 0 10000];

        R=20^2;
        K=P*H'/(H*P*H'+R);
        Ki=zeros(N,2);
        Pi=zeros(N,1);
        Piv=zeros(N,1);
        for i=1:N
            if isnan(Z(i))
                Xi=phi*Xi;
                P=phi*P*phi'+Q;
            else
                Xi=phi*Xi;
                P=phi*P*phi'+Q;
                Xi=Xi+K*(Z(i)-H*Xi);
                K=P*H'/(H*P*H'+R);
                Ki(i,:)=K';
                P=(eye(2)-K*H)*P;
            end

            xi(i)=Xi(1);
            
            XiE=phi*Xi;
            xiE1(i)=XiE(1);
            if i<N
                ErrE1(M,i)=(XiE(1)-X(i+1))^2;
                ErrVE1(M,i)=(XiE(2)-V(i+1))^2;
            end
            
            XiE=(phi^6)*XiE;
            xiE7(i)=XiE(1);
            if i<N-6
                ErrE7(M,i)=(XiE(1)-X(i+7))^2;
                ErrVE7(M,i)=(XiE(2)-V(i+7))^2;
            end
            

            Pi(i)=sqrt(P(1,1));
            Piv(i)=sqrt(P(2,2));
            Err(M,i)=(Xi(1)-X(i))^2;
            ErrV(M,i)=(Xi(2)-V(i))^2;
        end
    end

    ErrAvg=zeros(1,N-2);
    for j=1:N-2;
        ErrAvg(j)=sqrt((1/(M-1))*sum(Err(:,j+2)));
    end

    ErrAvgV=zeros(1,N-2);
    for j=1:N-2;
        ErrAvgV(j)=sqrt((1/(M-1))*sum(ErrV(:,j+2)));
    end
    
    ErrAvgE7=zeros(1,N-7);
    for j=1:N-7;
        ErrAvgE7(j)=sqrt((1/(M-1))*sum(ErrE7(:,j)));
    end
    
    ErrAvgE1=zeros(1,N-1);
    for j=1:N-1;
        ErrAvgE1(j)=sqrt((1/(M-1))*sum(ErrE1(:,j)));
    end
    
    ErrAvgVE7=zeros(1,N-7);
    for j=1:N-7;
        ErrAvgVE7(j)=sqrt((1/(M-1))*sum(ErrVE7(:,j)));
    end
    
    ErrAvgVE1=zeros(1,N-1);
    for j=1:N-1;
        ErrAvgVE1(j)=sqrt((1/(M-1))*sum(ErrVE1(:,j)));
    end
    
    plot(X)
    hold on
    plot(Z,'r')
    plot(xi,'g')
    plot(2:N+1,xiE1,'m')
    plot(8:N+7,xiE7,'k')
    title('Last run Trajectories plot')
    legend({'True trajectory','measurements','Kalman estimation', ...
        '1-step ahead extrapolation','7-steps ahead extrapolation'} ...
        ,'location','northeastoutside')
    grid on
    xlabel('Points')
    ylabel('Value')
    set(gcf,'position',[0,0,800,500]);
    
    figure
    subplot(1,2,1)
    plot(3:N,ErrAvg)
    hold on
    plot(Pi,'r')
    title('True estimation error Vs. calculation error')
    legend('True estimation error','Calculation error (Standard deviation)')
    xlabel('Points')
    ylabel('Value')
    grid on
    subplot(1,2,2)
    plot(2:N,ErrAvgE1,'r')
    hold on
    plot(8:N,ErrAvgE7)
    title('Extrapolation True estimation error')
    legend('1-step ahead extrapolation','7-step Ahead Extrapolation')
    xlabel('Points')
    ylabel('Value')
    grid on
    suptitle('X-position errors')
    set(gcf,'position',[0,0,800,500]);
    
    figure
    subplot(1,2,1)
    plot(3:N,ErrAvgV)
    hold on
    plot(Piv,'r')
    title('True estimation Vs. calculation error')
    legend('True estimation error','Calculation error (Standard deviation)')
    xlabel('Points')
    ylabel('Value')
    grid on
    subplot(1,2,2)
    plot(2:N,ErrAvgVE1,'r')
    hold on
    plot(8:N,ErrAvgVE7)
    title('Extrapolation true estimation error')
    legend('1-step ahead extrapolation','7-step Ahead Extrapolation')
    xlabel('Points')
    ylabel('Value')
    suptitle('Velocity errors')
    set(gcf,'position',[0,0,800,500]);
    grid on
end