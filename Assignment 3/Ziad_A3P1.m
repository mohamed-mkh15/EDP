close all;
clear;
normaldist=makedist('Normal',0,28);
W=random(normaldist,300,1);
X=zeros(300,1);
X(1)=10;
N=length(X);
for i=2:N
    X(i)=X(i-1)+W(i);
end

normaldist2=makedist('Normal',0,97);
eta=random(normaldist2,300,1);
Z=zeros(300,1);
N=length(X);
for i=1:N
    Z(i)=X(i)+eta(i);
end
plot(X)
hold on
plot(Z,'r')

sigmaw=28^2;
sigma_eta=97^2;

x=sigmaw/sigma_eta;
alpha=0.5*(-x+sqrt(x^2+4*x));

Xs=X;
for i=2:N
    Xs(i)=Xs(i-1)+alpha*(Z(i)-Xs(i-1));
end
figure
plot(Xs)
hold on
plot(X,'c')

Xsb=Xs;

for i=N-1:-1:1
   Xsb(i)=Xsb(i+1)+alpha*(Xs(i)-Xsb(i+1));
end

figure
plot(Xs)
hold on
plot(X,'c')
hold on
plot(Xsb,'r')