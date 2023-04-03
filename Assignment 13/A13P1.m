close all;
clear;

m=1;    %Extrapolation steps
N=500;  %Number of points
NE=N-m; %Number of points for the extrapolation error
Mm=500; %Number of runs

%Kalman errors initialization
ErrX=zeros(Mm,N);          %X true estimation error
ErrXE=zeros(Mm,NE);        %X-extrapolated true estimation error
ErrY=zeros(Mm,N);       %Y true estimation error
ErrYE=zeros(Mm,NE);     %Y-extrapolated true estimation error

for M=1:Mm
    
    load('theta.mat')
    
    X=zeros(1,N);           %X-position
    Y=zeros(1,N);           %Y-position
    X(1)=0;     %X initial value
    Y(1)=0;     %Y initial value
    V=10;       %velocity magnitude
    VX=zeros(1,N);          %X-Velocity
    VY=zeros(1,N);          %Y-Velocity
    VX(1)=V*cos(theta(1));     %X initial velocity
    VY(1)=V*sin(theta(1));     %Y initial velocity
    T=0.05;        %Time step
    
    sigma_a=1;      %Acceleration noise standard deviation
    
    normaldist=makedist('Normal',0,sigma_a);
    ax=random(normaldist,N,1);            %X-Acceleration noise
    
    normaldist=makedist('Normal',0,sigma_a);
    ay=random(normaldist,N,1);            %Y-Acceleration noise

    %Vectors generation
    for i=2:N
        X(i)=X(i-1)+VX(i-1)*T+0.5*ax(i-1)*T^2;        %X vector generation
        VX(i)=V*cos(theta(i));                    %X Velocity generation
        Y(i)=Y(i-1)+VY(i-1)*T+0.5*ay(i-1)*T^2;        %Y vector generation
        VY(i)=V*sin(theta(i));                    %Y Velocity generation        
    end
    
        
    sigma_x=3;           %X measurements noise standard deviation
    sigma_y=3;        %Y measurements noise standard deviation

    normaldist=makedist('Normal',0,sigma_x);
    eta_x=random(normaldist,1,N);       %X measurements noise vector

    normaldist=makedist('Normal',0,sigma_y);
    eta_y=random(normaldist,1,N);       %Y measurements noise vector
    
    zx=X+eta_x;
    zy=Y+eta_y;

    Z=[zx;zy];                %measurments

    %Kalman filter parameters initialization
    Xi=[zx(2);(zx(2)-zx(1))/T;zy(2);(zy(2)-zy(1))/T];   %State vector
    P=(10^4)*eye(4);         %P matrix

    %state space matrices
    phi=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];

    G=[0.5*T^2 0;T 0;0 0.5*T^2;0 T];
    
    sigma_a2=5;
    
    Q=G*G'*sigma_a2^2;
    
    R=[sigma_x^2 0;0 sigma_y^2];
    H=[1 0 0 0;0 0 1 0];
    %initial kalman gain
    K=P*H'/(H*P*H'+R);
    
    xi=zeros(N,1);
    yi=zeros(N,1);
    xiE=zeros(N,1);
    yiE=zeros(N,1);
    px=zeros(N,1);
    py=zeros(N,1);
    
    %Kalman filter
    for i=3:N
        Xi=phi*Xi;
        P=phi*P*phi'+Q;
        
        Xi=Xi+K*(Z(:,i)-H*Xi);
        K=P*H'/(H*P*H'+R);
        P=(eye(4)-K*H)*P;
        XiE=Xi;

        for mm=m
            XiE=phi*XiE;     %Extrapolated state vector
        end
        
        xi(i)=Xi(1);
        yi(i)=Xi(3);
        xiE(i)=XiE(1);
        yiE(i)=XiE(3);
        px(i)=P(1,1);
        py(i)=P(3,3);

        %True estimation error calculation
        ErrX(M,i)=(Xi(1)-X(i))^2;
        ErrY(M,i)=(Xi(3)-Y(i))^2;

        %Extrapolation error calculation
        if i<(N+m-1)
            ErrXE(M,i)=(XiE(1)-X(i+m))^2;
            ErrYE(M,i)=(XiE(3)-Y(i+m))^2;
        end
    end
end


%Trajectory plotting for visualization
figure(1)
plot(X,Y)
hold on
plot(xi,yi,'g')
plot(xiE,yiE,'m')
title('Trajectory polt')
legend('True trajectory','Filteration estimated','Extrapolation')

%plotting errors
figure(2)
subplot(1,2,1)
plotErr(ErrX,ErrXE,px,'X')
subplot(1,2,2)
plotErr(ErrY,ErrYE,py,'Y')
suptitle('True Estimation Error vs. Extrapolation error vs. standard deviation')