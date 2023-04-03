close all;
clear;

m=1;    %Extrapolation steps
N=500;  %Number of points
NE=N-m; %Number of points for the extrapolation error
Mm=500; %Number of runs

%Kalman errors initialization
ErrD=zeros(Mm,N);          %D true estimation error
ErrDE=zeros(Mm,NE);        %D-extrapolated true estimation error
Errbeta=zeros(Mm,N);       %beta true estimation error
ErrbetaE=zeros(Mm,NE);     %beta-extrapolated true estimation error

for M=1:Mm

    X=zeros(1,N);           %X-position
    Y=zeros(1,N);           %Y-position
    X(1)=1000;     %X initial value
    Y(1)=1000;     %Y initial value
    VX=zeros(1,N);          %X-Velocity
    VY=zeros(1,N);          %Y-Velocity
    VX(1)=100;     %X initial velocity
    VY(1)=100;     %Y initial velocity
    T=2;        %Time step
    
    sigma_a=0.3;      %Acceleration noise standard deviation
    
    normaldist=makedist('Normal',0,sigma_a);
    ax=random(normaldist,N,1);            %X-Acceleration noise
    
    normaldist=makedist('Normal',0,sigma_a);
    ay=random(normaldist,N,1);            %Y-Acceleration noise
    
    
    sigma_d_ob1=50;               %Range measurements noise standard deviation
    sigma_beta_ob1=0.004;          %Azimuth measurements noise standard deviation

    D=zeros(1,N);             %Range
    beta=zeros(1,N);          %Azimuth
    D(1)=sqrt(X(1)^2+Y(1)^2); %Range initial value
    beta(1)=atan2(X(1),Y(1));  %Azimuth initial value

    normaldist=makedist('Normal',0,sigma_d_ob1);
    eta_d_ob1=random(normaldist,N,1);            %Range measurements noise vector

    normaldist=makedist('Normal',0,sigma_beta_ob1);
    eta_b_ob1=random(normaldist,N,1);            %Azimuth measurements noise vector
    
    sigma_beta_ob2=0.001;          %Azimuth measurements noise standard deviation

    normaldist=makedist('Normal',0,sigma_beta_ob2);
    eta_b_ob2=random(normaldist,N,1);            %Azimuth measurements noise vector
    
    D_m=D;                                   %Range measurements
    beta_m=beta;                             %Azimuth measurements

    D_m(1)=D(1)+eta_d_ob1(1);                    %Range measurements initialization
    beta_m(1)=beta(1)+eta_b_ob1(1);              %Azimuth measurements initialization

    %Vectors generation
    for i=2:N
        X(i)=X(i-1)+VX(i-1)*T+0.5*ax(i-1)*T^2;        %X vector generation
        VX(i)=VX(i-1)+ax(i-1)*T;                    %X Velocity generation
        Y(i)=Y(i-1)+VY(i-1)*T+0.5*ay(i-1)*T^2;        %Y vector generation
        VY(i)=VY(i-1)+ay(i-1)*T;                    %Y Velocity generation        
        D(i)=sqrt(X(i)^2+Y(i)^2);     %Range vector generation
        beta(i)=atan2(X(i),Y(i));      %Azimuth vector generation
        
        if rem(i,2)==0
            D_m(i)=NaN;         %Range measurements vector generation
            beta_m(i)=beta(i)+eta_b_ob2(i);   %Azimuth measurements vector generation
        else
            D_m(i)=D(i)+eta_d_ob1(i);         %Range measurements vector generation
            beta_m(i)=beta(i)+eta_b_ob1(i);   %Azimuth measurements vector generation
        end
    end

    Z=[D_m;beta_m];                %Cartesian measurments
    
    x1=D_m(1)*sin(beta_m(1));
    y1=D_m(1)*cos(beta_m(1));
    x3=D_m(3)*sin(beta_m(3));
    y3=D_m(3)*cos(beta_m(3));
    %Kalman filter parameters initialization
    Xi=[x3;(x3-x1)/(2*T);y3;(y3-y1)/(2*T)];   %State vector
    P=(10^4)*eye(4);         %P matrix

    %state space matrices
    phi=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];

    G=[0.5*T^2 0;T 0;0 0.5*T^2;0 T];
    
    Q=G*G'*sigma_a^2;
    
    xi=Xi(1);
    yi=Xi(3);
    x2y2=xi^2+yi^2;
    h_ob1=[sqrt(xi^2+yi^2);atan2(xi,yi)];
    hprime_ob1=[xi/sqrt(x2y2) 0 yi/sqrt(x2y2) 0; yi/x2y2 0 -xi/x2y2 0];
    
    h_ob2=atan2(xi,yi);
    hprime_ob2=[yi/x2y2 0 -xi/x2y2 0];

    R_ob1=[sigma_d_ob1^2 0;0 sigma_beta_ob1^2];     %Error covariance matrix
    R_ob2=sigma_beta_ob2^2;
    
    %initial kalman gain
    K_ob1=P*hprime_ob1'/(hprime_ob1*P*hprime_ob1'+R_ob1);
    K_ob2=P*hprime_ob2'/(hprime_ob2*P*hprime_ob2'+R_ob2);

    Di=zeros(N,1);       %filtered Range
    betai=zeros(N,1);    %filtered Azimuth

    DiE=zeros(N,1);     %Extrapolated Range
    betaiE=zeros(N,1);  %Extrapolated Azimuth

    Ki=zeros(1,N);      %array of K(1,1)

    %Kalman filter
    for i=4:N
        if rem(i,2)~=0
            Xi=phi*Xi;
            P=phi*P*phi'+Q;

            xi=Xi(1);
            yi=Xi(3);
            x2y2=xi^2+yi^2;
            h=[sqrt(xi^2+yi^2);atan2(xi,yi)];
            hprime=[xi/sqrt(x2y2) 0 yi/sqrt(x2y2) 0; yi/x2y2 0 -xi/x2y2 0]; 

            Xi=Xi+K_ob1*(Z(:,i)-h);
            R=R_ob1;
            K_ob1=P*hprime'/(hprime*P*hprime'+R);
            P=(eye(4)-K_ob1*hprime)*P;
        else
            Xi=phi*Xi;
            P=phi*P*phi'+Q;

            xi=Xi(1);
            yi=Xi(3);
            x2y2=xi^2+yi^2;
            h=atan2(xi,yi);
            hprime=[yi/x2y2 0 -xi/x2y2 0]; 

            Xi=Xi+K_ob2*(Z(2,i)-h);
            R=R_ob2;
            K_ob2=P*hprime'/(hprime*P*hprime'+R);
            P=(eye(4)-K_ob2*hprime)*P;
        end
        XiE=Xi;

            for mm=m
                XiE=phi*XiE;     %Extrapolated state vector
            end


            %Range and azimuth estimation (filtered and extrapolated)
            Di(i)=sqrt(Xi(1)^2+Xi(3)^2);
            betai(i)=atan2(Xi(1),Xi(3));
            DiE(i)=sqrt(XiE(1)^2+XiE(3)^2);
            betaiE(i)=atan2(XiE(1),XiE(3));


            %True estimation error calculation
            ErrD(M,i)=(Di(i)-D(i))^2;
            Errbeta(M,i)=(betai(i)-beta(i))^2;

            %Extrapolation error calculation
            if i<N
                ErrDE(M,i)=(DiE(i)-D(i+m))^2;
                ErrbetaE(M,i)=(betaiE(i)-beta(i+m))^2;
            end
    end
end

%Trajectory plotting for visualization
polar(beta,D)
hold on
I = ~isnan(beta_m) & ~isnan(D_m);
polar(beta_m(I),D_m(I),'r')
polar(betai,Di,'g')
polar(betaiE,DiE,'m')
title('Trajectory polt using Range and azimuth')
legend({'True trajectory','Measurements (Observer 1)','Filteration estimated','Extrapolation'},...
    'location','southwest')

%plotting errors
figure
subplot(1,2,1)
plotErr(ErrD,ErrDE,sigma_d_ob1*ones(1,N),'Range')
subplot(1,2,2)
plotErr(Errbeta,ErrbetaE,sigma_beta_ob1*ones(1,N),'Azimuth')
plot(sigma_beta_ob2*ones(1,N),'g','linewidth',1.5)
suptitle('True Estimation Error vs. Extrapolation error vs. standard deviation')