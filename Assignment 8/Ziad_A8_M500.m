close all;
clear;
Err=zeros(500,200);
ErrV=zeros(500,200);
ErrS=zeros(500,200);
ErrSV=zeros(500,200);
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
    sigma_a=0.2^2;
    Q=G*G'*sigma_a;

    xi=zeros(1,200);
    Xi=[2;0];
    P=[10000 0; 0 10000];

    R=20^2;
    K=P*H'/(H*P*H'+R);
    Ki=zeros(N,2);
    Pi=zeros(N,1);
    Piv=zeros(N,1);
    xx=zeros(2,6);
    xiS=zeros(2,N);
    Parray1=cell(1,N);
    Parray2=cell(1,N);
    for i=1:N
        Xi=phi*Xi;
        P=phi*P*phi'+Q;
        Parray2{i}=P;
        Xi=Xi+K*(Z(i)-H*Xi);
        K=P*H'/(H*P*H'+R);
        Ki(i,:)=K';
        P=(eye(2)-K*H)*P;
        Parray1{i}=P;
        xi(i)=Xi(1);
        xiS(:,i)=Xi;
        Pi(i)=sqrt(P(1,1));
        Piv(i)=sqrt(P(2,2));
        Err(M,i)=(Xi(1)-X(i))^2;
        ErrV(M,i)=(Xi(2)-V(i))^2;
    end
    
    P=Parray1{end-1};
    PiSx=zeros(N-1,1);
    PiSv=zeros(N-1,1);
    
    for i=N-1:-1:1
        P1=Parray1{i};
        P2=Parray2{i+1};
        A=P1*phi'/P2;
        xiS(:,i)=xiS(:,i)+A*(xiS(:,i+1)-phi*xiS(:,i));
        P=P1+A*(P-P2)*A';
        PiSx(i)=sqrt(P(1,1));
        PiSv(i)=sqrt(P(2,2));
        ErrS(M,i)=(xiS(1,i)-X(i))^2;
        ErrSV(M,i)=(xiS(2,i)-V(i))^2;
    end
end

ErrAvg=zeros(1,N);
for j=1:N;
    ErrAvg(j)=sqrt((1/(M-1))*sum(Err(:,j)));
end

ErrAvgV=zeros(1,N);
for j=1:N-1;
    ErrAvgV(j)=sqrt((1/(M-1))*sum(ErrV(:,j)));
end

ErrAvgS=zeros(1,N-1);
for j=1:N-1;
    ErrAvgS(j)=sqrt((1/(M-1))*sum(ErrS(:,j)));
end

ErrAvgSV=zeros(1,N-1);
for j=1:N-1;
    ErrAvgSV(j)=sqrt((1/(M-1))*sum(ErrSV(:,j)));
end

subplot(2,2,1)
plot(ErrAvg)
hold on
plot(Pi,'r')
title('X')

subplot(2,2,2)
plot(1:N-1,ErrAvgS)
hold on
plot(PiSx,'r')
title('X smoothed')

subplot(2,2,3)
plot(ErrAvgV)
hold on
plot(Piv,'r')
title('V')

subplot(2,2,4)
plot(1:N-1,ErrAvgSV)
hold on
plot(PiSv,'r')
title('V smoothed')
