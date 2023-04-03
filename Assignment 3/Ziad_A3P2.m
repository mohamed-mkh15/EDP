close all;
clear;
normaldist=makedist('Normal',0,sqrt(10));
a=random(normaldist,300,1);
X=zeros(300,1);
X(1)=5;
V=zeros(300,1);
T=0.1;
N=length(X);
for i=2:N
    X(i)=X(i-1)+V(i-1)*T+a(i-1)*T^2/2;
    V(i)=V(i-1)++a(i-1)*T;
end

normaldist2=makedist('Normal',0,sqrt(500));
eta=random(normaldist2,300,1);
Z=zeros(300,1);
N=length(X);
for i=1:N
    Z(i)=X(i)+eta(i);
end

plot(X)
hold on
plot(Z,'r')

figure
k=1;
subplot(5,5,1)
plot(X)
hold on
plot(Z,'r')
for M=3:2:49
    k=k+1;
    XsR=X;
    j=(M-1)/2;
    for i=j+1:N-j
        XsR(i)=(1/M)*sum(Z((i-j):(i+j)));
    end
    XsR(1:j)=mean(X(1:j));
    XsR(end-j+1:end)=mean(X(end-j+1:end));
    subplot(5,5,k)
    plot(XsR)
    title(M);
end
k=1;
figure
subplot(4,4,1)
plot(X)
hold on
plot(Z,'r')
for alpha=0.01:0.01:0.15
    k=k+1;
    Xs=X;
    for i=2:N
        Xs(i)=Xs(i-1)+alpha*(Z(i)-Xs(i-1));
    end
    subplot(4,4,k)
    plot(Xs)
    title(alpha);
end
