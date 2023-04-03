function [] = modified_Kalman500(lambda1,string )
    %creating arrays to hold the values of the errors. 
    Err_P_x=ones(500,200); Err_P_v=ones(500,200); Err_P_a=ones(500,200);
    Err_F_x=ones(500,200); Err_F_v=ones(500,200); Err_F_a=ones(500,200);
    Err_Ext_x=ones(500,200); Err_Ext_v=ones(500,200); Err_Ext_a=ones(500,200);
    final_error=ones(3,198);
    final_error_p=ones(3,198);
    final_error_E=ones(3,198);
    for M=1:500
        T=1; 
        var_a=0.2^2;
        var_eta=20^2;
        var_zeta=var_a*(1-exp(-lambda1*T));
        N=200; 
        x=ones(1,N); x(1)=5; 
        a=ones(1,N); a(1)=normrnd(0,sqrt(var_zeta));
        z=ones(1,N); z(1)=x(1)+normrnd(0,sqrt(var_eta));
        v=ones(1,N); v(1)=0;
        for i=2:N
            a(i)=a(i-1)*exp(-lambda1*T)+normrnd(0,sqrt(var_zeta));
            x(i)=x(i-1)+v(i-1)*T+(a(i-1)*T^2)/2; 
            v(i)=v(i-1)+a(i-1)*T; 
            z(i)=x(i)+normrnd(0,sqrt(var_eta));
        end
        phi=[1 T (T^2)/2;0 1 T; 0 0 exp(-lambda1*T)]; 
        G=[0; 0;1]; 
        H=[1 0 0] ; 
        X_init=[2;0;0];                     % intial value of the state vector X=[x; v].
        p_init=[10000 0 0;0 10000 0; 0 0 10000];         % intial value for p matrix. 
        Q=G*G'*var_a;                     % intial value for covariance Q. 
        R=20^2;                           % intial value for covariance R.
        K_init=p_init*H'/(H*p_init*H'+R); % initial value for the Kalman gain. 
        filtered=ones(1,N);               % Matrix for holding the Kalman filter results. 
        K=ones(3,N);                      % Matrix for hilding the galman gain. 
        m=2;                              % number of extrapolation steps. 
        extrapolation=ones(3,N);          % Extrapolation results. 
        for i=1:N
           X_pre=phi*X_init;                  % prdection of X 
           p_pre=phi*p_init*phi'+Q;           % predection of P
           X_imp=X_pre+K_init*(z(i)-H*X_pre); % Improving X with measurments.
           X_ext_m=(phi^m)*X_imp;             % Extrapolation for m steps. 
           extrapolation(:,i)=X_ext_m;        % Extrapolation saving
           K(:,i)=K_init;                     % kalman Gain array 
           K_init=p_pre*H'/(H*p_pre*H'+R);    % Improving & updating the kalman gain. 
           filtered(i)=X_imp(1,1);            % trajectory saving 
           p_imp=(eye(3)-K_init*H)*p_pre;     % Improving P 
           Err_P_x(M,i)=sqrt(p_imp(1,1)); Err_P_v(M,i)=sqrt(p_imp(2,2)); Err_P_a(M,i)=sqrt(p_imp(3,3));    % standard deviation of error 
           Err_F_x(M,i)=(x(i)-X_imp(1,1))^2; Err_F_v(M,i)=(v(i)-X_imp(2,1))^2;  Err_F_a(M,i)=(a(i)-X_imp(3,1))^2;   % estimation error 
           %Err_Ext_x(M,i)=(x(i+1)- X_ex(1,1))^2;  Err_Ext_v(M,i)=(v(i+1)- X_ex(2,1))^2;  Err_Ext_a(M,i)=(a(i+1)- X_ex(3,1))^2; % extrapolation Error.
           p_init=p_imp;                      % updating P
           X_init=X_imp;                      % updating X
        end
        for i=1:N-m+1
        Err_Ext_x(M,i)=(x(i+m-1)-extrapolation(1,i))^2;
        Err_Ext_v(M,i)=(v(i+m-1)-extrapolation(2,i))^2;
        Err_Ext_a(M,i)=(a(i+m-1)-extrapolation(3,i))^2;
        end
        
    end
    for i =1:N-2
        final_error(1,i)=sqrt((1/499)*sum(Err_F_x(:,i+2)));
        final_error_p(1,i)=(1/499)*sum(Err_P_x(:,i+2));
        final_error(2,i)=sqrt((1/499)*sum(Err_F_v(:,i+2)));
        final_error_p(2,i)=(1/499)*sum(Err_P_v(:,i+2));        
        final_error(3,i)=sqrt((1/499)*sum(Err_F_a(:,i+2)));
        final_error_p(3,i)=(1/499)*sum(Err_P_a(:,i+2));
        
    end
    for i =1:N-m+1
       final_error_E(1,i)=(1/499)*sum(Err_Ext_x(:,i));
       final_error_E(2,i)=(1/499)*sum(Err_Ext_v(:,i));
       final_error_E(3,i)=(1/499)*sum(Err_Ext_a(:,i));
    end
    
figure()
subplot(2,1,1)
plot(1:198,final_error(1,:),1:198,final_error_p(1,:));
s='true vs calculated estimation error X';
s2=strcat(s,string);
title(s2,'color','r');
xlabel('trajectory points');
ylabel('Value');
legend('True Error','calculated error from p matrix');
grid on
% size(final_error_E)
subplot(2,1,2)
%plot(1:198,final_error_E(1,:));
plot(3:198,final_error_E(1,3:198));
s=' Extrapolation Error X';
s2=strcat(s,string);
title(s2,'color','r');
xlabel('trajectory points');
ylabel('Value');
ylim([0 300])
grid on
figure() 
subplot(2,1,1)
plot(1:198,final_error(2,:),1:198,final_error_p(2,:));
s='true vs calculated estimation error V';
s2=strcat(s,string);
title(s2,'color','r');
xlabel('trajectory points');
ylabel('Value');
legend('True Error','calculated error from p matrix');
grid on
subplot(2,1,2)
plot(3:198,final_error_E(2,3:198));
s=' Extrapolation Error V';
s2=strcat(s,string);
title(s2,'color','r');
%title("Figure(" + (graph_no+3) + ") Extrapolation Error V",'color','r');
xlabel('trajectory points');
ylabel('Value');
ylim([0 300])
grid on
figure()
subplot(2,1,1)
plot(1:198,final_error(3,:),1:198,final_error_p(3,:));
s='true vs calculated estimation error a';
s2=strcat(s,string);
title(s2,'color','r');
%title("Figure(" + (graph_no+4) + ") true vs calculated estimation error a",'color','r');
xlabel('trajectory points');
ylabel('Value');
legend('True Error','calculated error from p matrix');
grid on
subplot(2,1,2)
plot(3:198,final_error_E(3,3:198));
s='Extrapolation Error A';
s2=strcat(s,string);
title(s2,'color','r');
%title("Figure(" + (graph_no+5) + ") Extrapolation Error A",'color','r');
xlabel('trajectory points');
ylabel('Value');
ylim([0 300])
grid on
end

