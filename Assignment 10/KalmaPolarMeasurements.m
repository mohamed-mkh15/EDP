function [ Errx,ErrxE,Erry,ErryE,ErrD,ErrDE,Errbeta,ErrbetaE,Px,Py,...
    X,Y,xi,yi,xiE,yiE,X_m,Y_m ...
    ,D,beta,D_m,beta_m,Di,betai,DiE,betaiE,Ki,K500,CN] = KalmaPolarMeasurements( x1,y1,N,m,Mm,sigma_beta,sigma_d )

%   Detailed explanation goes here

    %Kalman errors initialization
    NE=N-m;
    Errx=zeros(Mm,N);          %X true estimation error
    ErrxE=zeros(Mm,NE);        %X-extrapolated true estimation error
    Erry=zeros(Mm,N);          %Y true estimation error
    ErryE=zeros(Mm,NE);        %Y-extrapolated true estimation error
    ErrD=zeros(Mm,N);          %D true estimation error
    ErrDE=zeros(Mm,NE);        %D-extrapolated true estimation error
    Errbeta=zeros(Mm,N);       %beta true estimation error
    ErrbetaE=zeros(Mm,NE);     %beta-extrapolated true estimation error
    K500=zeros(Mm,N);          %kalman filter gain matrix

    for M=1:Mm

        X=zeros(1,N);           %X-position
        Y=zeros(1,N);           %Y-position
        X(1)=x1;     %X initial value
        Y(1)=y1;     %Y initial value
        Vx=-50;     %X velocity
        Vy=-45;     %Y velocity
        T=2;        %Time step

        D=zeros(1,N);             %Range
        beta=zeros(1,N);          %Azimuth
        D(1)=sqrt(X(1)^2+Y(1)^2); %Range initial value
        beta(1)=atan2(X(1),Y(1));  %Azimuth initial value

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
            beta(i)=atan2(X(i),Y(i));      %Azimuth vector generation

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

        %Initial measurement error covariance matrix
        R=[(sin(beta_m(1))*sigma_d)^2+(D_m(1)*cos(beta_m(1))*sigma_beta)^2 ...
            sin(beta_m(1))*cos(beta_m(1))*(sigma_d^2-(D_m(1)*sigma_beta)^2); ...
            sin(beta_m(1))*cos(beta_m(1))*(sigma_d^2-(D_m(1)*sigma_beta)^2) ...
            (cos(beta_m(1))*sigma_d)^2+(D_m(1)*sin(beta_m(1))*sigma_beta)^2];
        K=P*H'/(H*P*H'+R);

        Px=zeros(N,1);       %X Calculation error initialization
        Py=zeros(N,1);       %Y Calculation error initialization
        xiE=zeros(N,1);      %Extrapolated X
        yiE=zeros(N,1);      %Extrapolated Y

        Di=zeros(N,1);       %filtered Range
        betai=zeros(N,1);    %filtered Azimuth

        DiE=zeros(N,1);     %Extrapolated Range
        betaiE=zeros(N,1);  %Extrapolated Azimuth

        Ki=zeros(1,N);      %array of K(1,1)
        
        CN=zeros(1,N);

        %Kalman filter
        for i=1:N

            Ki(i)=K(1,1);
            R=[(sin(beta_m(i))*sigma_d)^2+(D_m(i)*cos(beta_m(i))*sigma_beta)^2 ...
                sin(beta_m(i))*cos(beta_m(i))*(sigma_d^2-(D_m(i)*sigma_beta)^2); ...
                sin(beta_m(i))*cos(beta_m(i))*(sigma_d^2-(D_m(i)*sigma_beta)^2) ...
                (cos(beta_m(i))*sigma_d)^2+(D_m(i)*sin(beta_m(i))*sigma_beta)^2];
            CN(i)=cond(R);

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
            betai(i)=atan2(Xi(1),Xi(3));
            DiE(i)=sqrt(XiE(1)^2+XiE(3)^2);
            betaiE(i)=atan2(XiE(1),XiE(3));


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
        K500(M,:)=Ki;
    end
end

