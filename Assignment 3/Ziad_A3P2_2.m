close all;
clear;
normaldist=makedist('Normal',0,0.08);
w=random(normaldist,200,1);
X=zeros(200,1);
A=zeros(200,1);
T=32;
omega=2*pi/T;
A(1)=1;
X(1)=A(1)*sin(omega+3);
N=length(X);
for i=2:N
    A(i)=A(i-1)+w(i);
    X(i)=A(i)*sin(omega*i+3);
end

normaldist2=makedist('Normal',0,sqrt(0.05));
eta=random(normaldist2,200,1);
Z=zeros(200,1);

for i=1:N
    Z(i)=X(i)+eta(i);
end

plot(X)
hold on
plot(Z,'r')

%%
M=13;
j=(M-1)/2;
XsR=Z;
for i=j+1:N-j
    XsR(i)=(1/M)*sum(Z((i-j):(i+j)));
end
XsR(1:j)=mean(X(1:j));
XsR(end-j+1:end)=mean(X(end-j+1:end));

plot(XsR,'g')
