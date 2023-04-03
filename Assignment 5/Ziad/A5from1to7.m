close all;
clear;
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
plot(X)
hold on
plot(Z,'r')

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
xx=zeros(2,6);
xiE=zeros(N,1);
for i=1:N
    Xi=phi*Xi;
    Q=G*G'*sigma_a;
    P=phi*P*phi'+Q;
    Xi=Xi+K*(Z(i)-H*Xi);
    xi(i)=Xi(1);
    K=P*H'/(H*P*H'+R);
    Ki(i,:)=K';
    P=(eye(2)-K*H)*P;
    Pi(i)=sqrt(P(1,1));
    %Extrapolation
    XiE=Xi;
    for m=1:6
        XiE=phi*XiE;
        xx(:,m)=phi*XiE;
    end
    xiE(i)=XiE(1);
    figure(2)
    plot(i+(1:m),xx(1,:))
    plot(i,Xi(1),'go')
    hold on
end
figure()
plot(xi,'g')
plot(7:N+6,xiE,'k')
figure()
plot(xi,'r')
figure ()
plot(Ki(:,1))
hold on 
plot(Ki(:,2),'r')

figure()
plot(Pi)
figure()
ylimm=get(gca,'ylim');
xlimm=get(gca,'xlim');
figure()
set(gca,'ylim',ylimm);
set(gca,'xlim',xlimm);