close all;
clear;
Err=zeros(500,200);
ErrV=zeros(500,200);
Erra=zeros(500,200);
for M=1:500
    X=zeros(1,200);
    X(1)=5;
    V=zeros(1,200);
    V(1)=1;
    T=1;
    N=length(X);

    lambda=0.1;
    elambda=exp(-lambda*T);
    elambda1=elambda;

    sigma_a=0.2;
    sigma_zeta=sqrt((sigma_a^2)*(1-elambda^2));
    normaldist=makedist('Normal',0,sigma_zeta);
    zeta=random(normaldist,200,1);

    a=zeros(1,N);
    a(1)=normrnd(0,sigma_a);

    for i=2:N
        a(i)=elambda*a(i-1)+zeta(i);
    end

    for i=2:N
        X(i)=X(i-1)+V(i-1)*T+a(i-1)*T^2/2;
        V(i)=V(i-1)+a(i-1)*T;
    end
    
    normaldist2=makedist('Normal',0,20);
    eta=random(normaldist2,200,1);

    Z=zeros(1,200);

    for i=1:N
        Z(i)=X(i)+eta(i);
    end
    X1=X;

    phi=[1 T T^2/2;0 1 T;0 0 elambda1];
    G=[0; 0;1];
    H=[1 0 0];

    xi=zeros(1,200);
    Xi=[2;0;0];
    P=10000*eye(3);

    R=20^2;
    sigma_a=0.2^2;
    K=P*H'/(H*P*H'+R);
    Ki=zeros(N,3);
    Px=zeros(N,1);
    Pv=zeros(N,1);
    Pa=zeros(N,1);

    for i=1:2
        Xi=phi*Xi;
        Q=G*G'*sigma_a;
        P=phi*P*phi'+Q;
        Xi=Xi+K*(Z(i)-H*Xi);
        xi(i)=Xi(1);
        K=P*H'/(H*P*H'+R);
        Ki(i,:)=K';
        P=(eye(3)-K*H)*P;
        Px(i)=sqrt(P(1,1));
        Pv(i)=sqrt(P(2,2));
        Pa(i)=sqrt(P(3,3));
    end
    for i=3:N
        Xi=phi*Xi;
        Q=G*G'*sigma_a;
        P=phi*P*phi'+Q;
        Xi=Xi+K*(Z(i)-H*Xi);
        xi(i)=Xi(1);
        K=P*H'/(H*P*H'+R);
        Ki(i,:)=K';
        P=(eye(3)-K*H)*P;
        Px(i)=sqrt(P(1,1));
        Pv(i)=sqrt(P(2,2));
        Pa(i)=sqrt(P(3,3));
        Err(M,i)=(Xi(1)-X(i))^2;
        ErrV(M,i)=(Xi(2)-V(i))^2;
        Erra(M,i)=(Xi(3)-a(i))^2;
     end
end
ErrAvg=zeros(1,N-2);
for j=1:N-2;
    ErrAvg(j)=sqrt((1/(M-1))*sum(Err(:,j+2)));
end
figure(1)
plot(3:N,ErrAvg)
hold on
plot(Px,'r')
title('error x')

figure(2)
ErrAvgV=zeros(1,N-2);
for j=1:N-2;
    ErrAvgV(j)=sqrt((1/(M-1))*sum(ErrV(:,j+2)));
end
plot(3:N,ErrAvgV)
hold on
plot(Pv,'r')
title('error v')

figure(3)
ErrAvga=zeros(1,N-2);
for j=1:N-2;
    ErrAvga(j)=sqrt((1/(M-1))*sum(Erra(:,j+2)));
end
plot(3:N,ErrAvga)
hold on
plot(Pa,'r')
error(a)

figure(4)
plot(X)
hold on
plot(xi,'g')
plot(Z,'r')
