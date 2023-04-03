function [] = Kalman500(lambda1,lambda2)
    Err_P=ones(500,200);
    Err_F=ones(500,200);
    final_error=ones(1,198);
    final_error_p=ones(1,198);
    for M=1:500
        T=1; 
        var_a=0.2^2;
        var_eta=20^2;
        var_zeta=var_a*(1-exp(-lambda1*T));
        var_vsigma=var_eta*(1-exp(-lambda2*T));
        N=200; 
        x=ones(1,N); x(1)=5; 
        a=ones(1,N); a(1)=normrnd(0,sqrt(var_zeta));
        eta=ones(1,N); eta(1)=normrnd(0,sqrt(var_vsigma));
        z=ones(1,N); z(1)=x(1)+eta(1);
        v=ones(1,N); v(1)=0;
        for i=2:N
            a(i)=a(i-1)*exp(-lambda1*T)+normrnd(0,sqrt(var_zeta));
            eta(i)=eta(i-1)*exp(-lambda2*T)+normrnd(0,sqrt(var_vsigma));
            x(i)=x(i-1)+v(i-1)*T+(a(i-1)*T^2)/2; 
            v(i)=v(i-1)+a(i-1)*T; 
            z(i)=x(i)+eta(i);
        end
        phi=[1 T; 0 1]; 
        G=[(T^2)/2; T]; 
        H=[1 0] ; 
        X_init=[2;0];                     % intial value of the state vector X=[x; v].
        p_init=[10000 0;0 10000];         % intial value for p matrix. 
        Q=G*G'*var_a;                     % intial value for covariance Q. 
        R=20^2;                           % intial value for covariance R.
        K_init=p_init*H'/(H*p_init*H'+R); % initial value for the Kalman gain. 
        filtered=ones(1,N);               % Matrix for holding the Kalman filter results. 
        K=ones(2,N);                      % Matrix for hilding the galman gain. 
        Err_p=ones(1,N);                  % Aray for holding the error standardd deviation.  
        Err_f=ones(1,N);                  % estimation error
        for i=1:N
           X_pre=phi*X_init;                  % prdection of X 
           p_pre=phi*p_init*phi'+Q;           % predection of P
           X_imp=X_pre+K_init*(z(i)-H*X_pre); % Improving X with measurments. 
           K(:,i)=K_init;                     % kalman Gain array 
           K_init=p_pre*H'/(H*p_pre*H'+R);    % Improving & updating the kalman gain. 
           filtered(i)=X_imp(1,1);            % trajectory saving 
           p_imp=(eye(2)-K_init*H)*p_pre;     % Improving P 
           Err_p(i)=sqrt(p_imp(1,1));         % standard deviation of error 
           Err_f(i)=(x(i)-X_imp(1,1))^2;      % estimation error 
           p_init=p_imp;                      % updating P
           X_init=X_imp;                      % updating X
        end
        Err_P(M,:)=Err_p;
        Err_F(M,:)=Err_f;
    end
    for i =1:N-2
        final_error(i)=sqrt((1/499)*sum(Err_F(:,i+2)));
        final_error_p(i)=(1/499)*sum(Err_P(:,i+2));
    end

figure()
plot(1:198,final_error,1:198,final_error_p);
title('true vs standard deviation of estimation error','color','r');
xlabel('trajectory points');
ylabel('Value');
legend('True Error','calculated error from p matrix');
grid on
end

