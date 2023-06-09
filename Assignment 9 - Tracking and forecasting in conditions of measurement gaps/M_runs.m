function [] = M_runs(no_of_runs,min_probability)
    N=200;              % size of trajectory
    err_run_X = zeros(N-2,no_of_runs);
    err_run_V = zeros(N-2,no_of_runs);
    err_run_X_ex = zeros(N-2,no_of_runs);
    err_run_V_ex = zeros(N-2,no_of_runs);
    
    for r = 1:no_of_runs
        x = zeros(N,1);     % real data for position
        V = zeros(N,1);     % real data for velocity
        a = zeros(N,1);     % noise
        z = zeros(N,1);     % measurements
        s_alpha = 0.2^2;    % sigma_w^2
        s_eta = 20^2;

        x(1) = 5;           % initial position
        V(1) = 1;           % initial velocity
        a(1) = normrnd(0,sqrt(s_eta));
        z(1) = x(1) + a(1);
        T = 1;              % time step

        for i = 2:N
            a(i) = normrnd(0,sqrt(s_alpha));
            x(i) = x(i-1) + V(i-1)*T + a(i)*T^2 /2; 
            V(i) = V(i-1) + a(i)*T;  
            z(i) = x(i) + normrnd(0,sqrt(s_eta));
        end


        phi = [1 T;0 1];    % transition matrix that relates Xi to Xi-1
        G = [T^2/2;T];      % input matrix, that determines how random acceleration affects state vector
        H = [1 0];          % observation matrix
        Q = G*G'*s_alpha;   % covariance matrix of state noise
        K = zeros(N,2);     %  array to save filter gain


        I = eye(2);
        R = s_eta;          % covariance matrix of measurements noise

        X = zeros(N,2);     % state vector, that describes full state of the system
        x_extrapolated = zeros(N-6,1);
        X(1,:)=[2,0];                    % Initial filtered estimate
        P_filteration = [10000 0; 0 10000];          % Initial filtration error covariance matrix
        P_value_X = zeros(N,1);            % array to save the caculation error in X
        P_value_V = zeros(N,1);            % array to save the caculation error in V
        K(1,:) = (P_filteration*H'*(H*P_filteration*H'+R)^-1)';  % initial value of filter gain

        m = 7;                  % extrapolation steps
        P_filteration_cell = cell(N,1);     % cell to save the P_filteration matrices
        P_prediction_cell = cell(N,1);     % cell to save the P_prediction matrices

        % kalman filter
        for i = 2:N
            % prediction(extrapolation)
                X_predicted = phi * X(i-1,:)';      % Prediction of state vector at time i using i ? 1 measurements
                P_prediction = phi*P_filteration*phi' + Q;                 % Prediction error covariance matrix
            % estimation
                X(i,:) = X_predicted + K(i-1,:)'*(z(i)- H*X_predicted);     % Improved estimate by incorporating a new measurement
                K_updated = P_prediction*H'*(H*P_prediction*H'+R)^-1;     % Filter gain, weight of residual
                P_filteration = (I - K_updated*H)*P_prediction;            % Filtration error covariance matrix
            % gaps
                probability = rand(1);
                if probability <= min_probability
                    X(i,:) = X_predicted;
                    P_filteration = P_prediction;
                end
            % extrapolation
            if i >= m
                x1 = phi^(m-1) * X(i-m+1,:)';
                x_extrapolated(i-m+1)=x1(1);
            end
            % save values to be used in the smoothing step
                P_filteration_cell{i} = P_filteration; 
                P_prediction_cell{i} = P_prediction; 

        P_value_X(i)=sqrt(P_filteration(1,1)); 
        P_value_V(i)=sqrt(P_filteration(2,2));        
        K(i,:) = K_updated;

        end
        
        err_run_X(:,r)=(x(3:N)- X(3:N,1)).^2;
        err_run_V(:,r)=(V(3:N)- X(3:N,2)).^2;   
    
    end
    
    error_true_X = zeros(N-2,1);
    error_true_V = zeros(N-2,1);

    for i = 1:N-3
        error_true_X(i) = sqrt((1/(no_of_runs-1)*sum(err_run_X(i,:))));
        error_true_V(i) = sqrt((1/(no_of_runs-1)*sum(err_run_V(i,:))));
    end

        error_true_X(N-2) = sqrt((1/(no_of_runs-1)*sum(err_run_X(N-2,:))));
        error_true_V(N-2) = sqrt((1/(no_of_runs-1)*sum(err_run_V(N-2,:))));

    
    % Plot true & filteration error for X
    figure()
    subplot(2,1,1)
    plot(3:N,error_true_X); hold on
    plot(3:N,P_value_X(3:N),'r')
    title('True error vs filteration error in X','color','r')
    xlabel('iteration')
    legend('true error','filteration error')
    ylabel('value')       

    
     % Plot true & filteration error for V
    figure()
    %subplot(2,1,1)
    plot(3:N,error_true_V); hold on
    plot(3:N,P_value_V(3:N),'r')
    title('True error vs filteration error in V','color','r')
    xlabel('iteration')
    legend('true error','filteration error')
    ylabel('value') 
    
  
    
end