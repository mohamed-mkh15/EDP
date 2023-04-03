close all;
clear;
Err=zeros(500,200);
ErrE=zeros(500,194);
for M=1:500
    X=zeros(1,200);
    X(1)=5;
    V=zeros(1,200);
    V(1)=1;
    T=1;
    N=length(X);

    lambda=0.1;
    elambda=exp(-lambda*T);

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
    
    lambda=1000;
    elambda=exp(-lambda*T);

    sigma_eta=20;
    sigma_zeta2=sqrt((sigma_eta^2)*(1-elambda^2));
    normaldist=makedist('Normal',0,sigma_zeta2);
    zeta2=random(normaldist,200,1);

    eta=zeros(1,N);
    eta(1)=normrnd(0,sigma_eta);

    for i=2:N
        eta(i)=elambda*eta(i-1)+zeta2(i);
    end

    Z=zeros(1,200);

    for i=1:N
        Z(i)=X(i)+eta(i);
    end
    X1=X;

    phi=[1 T;0 1];
    G=[T^2; T];
    H=[1 0];

    xi=zeros(1,200);
    Xi=[2;0];
    P=[10000 0; 0 10000];

    R=20^2;
    sigma_a=0.2^2;
    K=P*H'/(H*P*H'+R);
    Ki=zeros(N,2);
    Pi=zeros(N,1);

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
     end
end
ErrAvg=zeros(1,N-2);
for j=1:N-2;
    ErrAvg(j)=sqrt((1/(M-1))*sum(Err(:,j+2)));
end
plot(3:N,ErrAvg)
hold on
plot(Pi,'r')

figure
plot(X)
hold on
plot(xi,'g')
plot(Z,'r')
