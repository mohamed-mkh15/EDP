close all;
clear;

m=1;    %Extrapolation steps
N=26;   %Number of points
NE=N-m; %Number of points for the extrapolation error

Errx=zeros(500,N);          %X true estimation error
ErrxE=zeros(500,NE);        %X-extrapolated true estimation error
Erry=zeros(500,N);          %Y true estimation error
ErryE=zeros(500,NE);        %Y-extrapolated true estimation error
ErrD=zeros(500,N);          %D true estimation error
ErrDE=zeros(500,NE);        %D-extrapolated true estimation error
Errbeta=zeros(500,N);       %beta true estimation error
ErrbetaE=zeros(500,NE);     %beta-extrapolated true estimation error

for M=1:500

    X=zeros(1,N);           %X-position
    Y=zeros(1,N);           %Y-position
    X(1)=13500/sqrt(2);     %X initial value
    Y(1)=13500/sqrt(2);     %Y initial value
    Vx=-50;     %X velocity
    Vy=-45;     %Y velocity
    T=2;        %Time step

    D=zeros(1,N);             %Range
    beta=zeros(1,N);          %Azimuth
    D(1)=sqrt(X(1)^2+Y(1)^2); %Range initial value
    beta(1)=atan(X(1)/Y(1));  %Azimuth initial value

    sigma_d=20;               %Range measurements noise standard deviation
    sigma_beta=0.02;          %Azimuth measurements noise standard deviation
    
    normaldist=makedist('Normal',0,sigma_d);
    eta_d=random(normaldist,N,1);            %Range measurements noise vector

    normaldist=makedist('Normal',0,sigma_beta);
    eta_b=random(normaldist,N,1);            %Azimuth measurements noise vector
    D_m=D;                                   %Range measurements
    beta_m=beta;                             %Azimuth measurements

    D_m(1)=D(1)+eta_d(1);                    %Range measurements initialization
    beta_m(1)=beta(1)+eta_b(1);              %Azimuth measurements initialization
    
    %Vectors generation
    for i=2:N
        X(i)=X(i-1)+Vx*T;             %X vector generation
        Y(i)=Y(i-1)+Vy*T;             %Y vector generation
        D(i)=sqrt(X(i)^2+Y(i)^2);     %Range vector generation
        beta(i)=atan(X(i)/Y(i));      %Azimuth vector generation

        D_m(i)=D(i)+eta_d(i);         %Range measurements vector generation
        beta_m(i)=beta(i)+eta_b(i);   %Azimuth measurements vector generation
    end

    X_m=D_m.*sin(beta_m);       %X measurements vector generation    
    Y_m=D_m.*cos(beta_m);       %Y measurements vector generation

    Z=[X_m;Y_m];                %Cartesian measurments

    %Kalman filter parameters initialization
    xi=zeros(1,N);
    yi=zeros(1,N);
    Xi=[40000;-20;40000;-20];   %State vector
    P=(10^10)*eye(4);           %P matrix
    
    %state space matrices
    phi=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];
    H=[1 0 0 0;0 0 1 0];
    
    %measurement error covariance matrix
    R=[(sin(beta_m(1))*sigma_d)^2+(D_m(1)*cos(beta_m(1))*sigma_beta)^2 ...
        sin(beta_m(1))*cos(beta_m(1))*(sigma_d^2-(D_m(1)*sigma_beta)^2); ...
        sin(beta_m(1))*cos(beta_m(1))*(sigma_d^2-(D_m(1)*sigma_beta)^2) ...
        (cos(beta_m(1))*sigma_d)^2+(D_m(1)*sin(beta_m(1))*sigma_beta)^2];      

    K=P*H'/(H*P*H'+R);   %Kalman gain initialization

    Px=zeros(N,1);       %X Calculation error initialization
    Py=zeros(N,1);       %Y Calculation error initialization
    xiE=zeros(N,1);      %Extrapolated X
    yiE=zeros(N,1);      %Extrapolated Y

    Di=zeros(N,1);       %filtered Range
    betai=zeros(N,1);    %filtered Azimuth

    DiE=zeros(N,1);     %Extrapolated Range
    betaiE=zeros(N,1);  %Extrapolated Azimuth

    %Kalman filter
    for i=1:N

        R=[(sin(beta_m(i))*sigma_d)^2+(D_m(i)*cos(beta_m(i))*sigma_beta)^2 ...
            sin(beta_m(i))*cos(beta_m(i))*(sigma_d^2-(D_m(i)*sigma_beta)^2); ...
            sin(beta_m(i))*cos(beta_m(i))*(sigma_d^2-(D_m(i)*sigma_beta)^2) ...
            (cos(beta_m(i))*sigma_d)^2+(D_m(i)*sin(beta_m(i))*sigma_beta)^2];

        Xi=phi*Xi;
        P=phi*P*phi';
        Xi=Xi+K*(Z(:,i)-H*Xi);
        xi(i)=Xi(1);
        yi(i)=Xi(3);

        K=P*H'/(H*P*H'+R);
        P=(eye(4)-K*H)*P;
        XiE=Xi;
        
        for mm=m
            XiE=phi*XiE;     %Extrapolated state vector
        end
        
        %Extrapolated X,Y
        xiE(i)=XiE(1);
        yiE(i)=XiE(3);
        
        %Range and azimuth estimation (filtered and extrapolated)
        Di(i)=sqrt(Xi(1)^2+Xi(3)^2);
        betai(i)=atan(Xi(1)/Xi(3));
        DiE(i)=sqrt(XiE(1)^2+XiE(3)^2);
        betaiE(i)=atan(XiE(1)/XiE(3));
        
        %True estimation error calculation
        Errx(M,i)=(Xi(1)-X(i))^2;
        Erry(M,i)=(Xi(3)-Y(i))^2;
        ErrD(M,i)=(Di(i)-D(i))^2;
        Errbeta(M,i)=(betai(i)-beta(i))^2;
        
        %Covariance matrix error
        Px(i)=sqrt(P(1,1));
        Py(i)=sqrt(P(3,3));
        
        %Extrapolation error calculation
        if i<N
            ErrxE(M,i)=(XiE(1)-X(i+m))^2;
            ErryE(M,i)=(XiE(3)-Y(i+m))^2;
            ErrDE(M,i)=(DiE(i)-D(i+m))^2;
            ErrbetaE(M,i)=(betaiE(i)-beta(i+m))^2;
        end
    end
end

%Error plotting
subplot(2,2,1)
plotErr(Errx,ErrxE,Px,'X')
subplot(2,2,2)
plotErr(Erry,ErryE,Py,'Y')
subplot(2,2,3)
plotErr(ErrD,ErrDE,NaN,'D')
subplot(2,2,4)
plotErr(Errbeta,ErrbetaE,NaN,'Beta')

%Trajectory plotting for visualization
figure
subplot(1,2,1)
polar(beta,D)
hold on
polar(beta_m,D_m,'r')
polar(betai,Di,'g')
polar(betaiE,DiE,'m')
title('Trajectory polt using Range and azimuth')
legend('True trajectory','noise','Filteration estimated','Extrapolation')

subplot(1,2,2)
plot(X,Y)
hold on
plot(X_m,Y_m,'r')
plot(xi,yi,'g')
plot(xiE,yiE,'m')
title('Trajectory polt using Cartesian coordinates')
xlabel('X-position')
ylabel('Y-position')
legend('True trajectory','noise','Filteration estimated','Extrapolation')
grid on

betai_degree = zeros(length(xi),1);
for i = 1:length(xi)
    betai_degree(i) = rad2deg(betai(i));
end
figure
subplot(2,1,1)
plot(betai_degree,xi)
subplot(2,1,2)
plot(betai,xi)