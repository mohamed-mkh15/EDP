close all;
clear;
Err=zeros(500,200);
ErrE=zeros(500,194);
for M=1:500
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
    for i=1:N
        Z(i)=X(i)+eta(i);
    end
    X1=X;

    phi=[1 T;0 1];
    G=[T^2; T];
    H=[1 0];

    xi=zeros(1,200);
    Xi=[2;0];
    P=[100 0; 0 100];

    R=20^2;
    sigma_a=0.2^2;
    K=P*H'/(H*P*H'+R);
    Ki=zeros(N,2);
    Pi=zeros(N,1);
    xiE=zeros(N,1);
    xx=zeros(2,7);
    for i=1:2
        Xi=phi*Xi;
        Q=G*G'*sigma_a;
        P=phi*P*phi'+Q;
        Xi=Xi+K*(Z(i)-H*Xi);
        xi(i)=Xi(1);
        K=P*H'/(H*P*H'+R);
        Ki(i,:)=K';
        P=(eye(2)-K*H)*P;
        Pi(i)=sqrt(P(1,1));
        XiE=Xi;
        for m=1:6
            XiE=phi*XiE;
            xx(:,m)=phi*XiE;
        end
        xiE(i)=XiE(1);
        ErrE(M,i)=(XiE(1)-X(i+6))^2;
    end
     for i=3:N
        Xi=phi*Xi;
        Q=G*G'*sigma_a;
        P=phi*P*phi'+Q;
        Xi=Xi+K*(Z(i)-H*Xi);
        xi(i)=Xi(1);
        K=P*H'/(H*P*H'+R);
        Ki(i,:)=K';
        P=(eye(2)-K*H)*P;
        Pi(i)=sqrt(P(1,1));
        Err(M,i)=(Xi(1)-X(i))^2;
        XiE=Xi;
        for m=1:6
            XiE=phi*XiE;
            xx(:,m)=phi*XiE;
        end
        xiE(i)=XiE(1);
        if i<N-5
            ErrE(M,i)=(XiE(1)-X(i+6))^2;
        end
    end
end
ErrAvg=zeros(1,N-2);
for j=1:N-2;
    ErrAvg(j)=sqrt((1/(M-1))*sum(Err(:,j+2)));
end
ErrAvgE=zeros(1,N-6);
for j=1:N-6;
    ErrAvgE(j)=sqrt((1/(M-1))*sum(ErrE(:,j)));
end
plot(3:N,ErrAvg)
hold on
plot(Pi,'r')
figure
plot(7:N,ErrAvgE)
figure
plot(X)
hold on
plot(xi,'g')
plot(7:N+6,xiE,'k')