% Assignment 13 part 1
clc; clear 
load('theta.mat'); 
n = 500;        % size of trajectory
m = 1;          % no. of runs

% Crating Arrays to save the values of Errors over 500 runs. 
XF_Err = ones(500,n); XF_E = ones(1,n);
YF_Err = ones(500,n); YF_E = ones(1,n);
XP_Err = ones(500,n); XP_E = ones(1,n);
YP_Err = ones(500,n); YP_E = ones(1,n);
XE_Err = ones(500,n-m); XE_E = ones(1,n-m);
YE_Err = ones(500,n-m); YE_E = ones(1,n-m);

for j = 1:500
    N = length(theta);
    T = 0.05;
    x = zeros(1,N);
    y = zeros(1,N);
    sigma_z = 3;
    sigma_a = 5;
    v = 10;
    zx = zeros(1,N); zx(1) = normrnd(0,sigma_z);
    zy = zeros(1,N); zy(1) = normrnd(0,sigma_z);
    for i = 2:N
        x(i) = x(i-1)+v*cos(theta(i-1))*T+0.5*T*normrnd(0,1);
        y(i) = y(i-1)+v*sin(theta(i-1))*T+0.5*T*normrnd(0,1);
        zx(i) = x(i)+normrnd(0,sigma_z);
        zy(i) = y(i)+normrnd(0,sigma_z);
    end
    % plot(x,y,zx,zy,'*')
    XF = zeros(4,N);
    XE = zeros(4,N);
    X0 = [zx(2);
        (zx(2)-zx(1))/T;
        zy(2);
        (zy(2)-zy(1))/T];
    phi = [1 T 0 0;
        0 1 0 0;
        0 0 1 T;
        0 0 0 1];
    G = [(T^2)/2 0;
        T 0;
        0 (T^2)/2;
        0 T];
    H = [1 0 0 0;
        0 0 1 0];
    Q = G*G'.*sigma_a^2;
    P = diag([10^4 10^4 10^4 10^4 ]);
    R = [sigma_z^2 0;
        0 sigma_z^2];
    K = P*H'/(H*P*H'+R);
    for i = 3:N
        Z = [zx(i);
           zy(i)];
        Xp = phi*X0;
        Pex = phi*P*phi'+Q;
        Xf = Xp+K*(Z-H*Xp);
        Xe = phi*Xf;
        XF(:,i) = Xf;
        XE(:,i) = Xe;
        K = Pex*H'/(H*Pex*H'+R);
        P = (eye(4)-K*H)*Pex;
        X0 = Xf;
        % Errors Calculation. 
        XF_Err(j,i) = (x(i)-Xf(1))^2;
        YF_Err(j,i) = (y(i)-Xf(3))^2; 
        XP_Err(j,i) = sqrt(P(1,1)); 
        YP_Err(j,i) = sqrt(P(3,3));
        if i<N
            XE_Err(j,i) = (x(i+1)-Xe(1))^2;
            YE_Err(j,i) = (y(i+1)-Xe(3))^2;
        end
    end
    
end 
for i = 3:N
   XF_E(i) = sqrt((1/499)*sum(XF_Err(:,i))); 
   YF_E(i) = sqrt((1/499)*sum(YF_Err(:,i)));
   XP_E(i) = (1/499)*sum(XP_Err(:,i)); 
   YP_E(i) = (1/499)*sum(YP_Err(:,i)); 
   if i<N
       XE_E(i) = sqrt((1/499)*sum(XE_Err(:,i))); 
       YE_E(i) = sqrt((1/499)*sum(YE_Err(:,i))); 
   end
end 
plot(x,y,XF(1,:),XF(3,:),'g',XE(1,:),XE(3,:),'r',zx,zy,'*')
title('Trajectory for visualization','color','r')
xlabel('X-coordinate')
ylabel('y-coordinate')
legend('True trajectory','Filtered','Extapolated','measurment','location','northeastoutside')
figure()
plot(1:500-2,XF_E(3:end),1:500-2,XP_E(3:end));
hold on 
plot(2:500-2,XE_E(3:end));
hold off 
title('Errors of X coordinate estimation','color','r')
xlabel('iteration')
ylabel('value')
legend('True','Calculated','Extapolation','location','northeastoutside')
figure()
plot(1:500-2,YF_E(3:end),1:500-2,YP_E(3:end));
hold on 
plot(2:500-2,YE_E(3:end));
hold off 
title('Errors of Y coordinate estimation','color','r')
xlabel('iteration')
ylabel('value')
legend('True','Calculated','Extapolation','location','northeastoutside')

