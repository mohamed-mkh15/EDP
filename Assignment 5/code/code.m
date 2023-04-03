
clc; clear; close all;

N=200;              % size of the array
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
G = [T^2/2;T];      % input matrix, that determines how random acceleration ?? affects state vector
H = [1 0];          % observation matrix
Q = G*G'*s_alpha;
K = zeros(N,2);


I = eye(2);
R = s_eta;

X = zeros(N,2);          % state vector, that describes full state of the system
x_extrapolated = zeros(N-6,1);
X(1,:)=[2,0];       % Initial filtered estimate
P = [10000 0; 0 10000];   % Initial filtration error covariance matrix
P_vec = zeros(N,1);
K(1,:) = (P*H'*(H*P*H'+R)^-1)';

% kalman filter
for i = 2:N
    % prediction
        X_predicted = phi * X(i-1,:)';
        P = phi*P*phi' + Q;
    % estimation
        X(i,:) = X_predicted + K(i-1,:)'*(z(i)- H*X_predicted);
        K_updated = P*H'*(H*P*H'+R)^-1
        P = (I - K_updated*H)*P;
    % extrapolation
        if i >= 7
            x1 = phi^6 * X(i-6,:)'
            x_extrapolated(i-6)=x1(1);
        end
        
P_vec(i)=sqrt(P(1,1));        
K(i,:) = K_updated;

end

figure(1)
plot(x); hold on
plot(z); hold on
plot(X(:,1),'k')
legend('true','measurements','filterd data')
title ('Kalman Filter')

figure(2)
plot(K)
title('Filter gain K')

figure(3)
plot(P_vec)
title('Filteration error covariance P')

figure(4)
plot(X(:,1),'k'); hold on
plot(7:N,x_extrapolated)
title('Extrapolation on m=7 steps')
legend('filterred data','extrapolated data')



