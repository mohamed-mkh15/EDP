clear;
close all;
load('data_group5.mat');
Y=data(:,1);
M=data(:,2);
S=data(:,3);
T=Y+1/12*M;
plot(T,S)
hold on
SR=S;
N=length(S);
for i=7:N-6
    SR(i)=1/24*(S(i-6)+S(i+6))+1/12*(S(i-5)+S(i-4)+S(i-3)+S(i-2)+S(i-1)+S(i)+S(i+1)+S(i+2)+S(i+3)+S(i+4)+S(i+5));
end
SR(1:6)=mean(S(1:6));
SR(N-5:N)=mean(S(N-5:N));

plot(T,SR,'r')

[ Id_R,Iv_R ] = IdIv( S,SR );

%%
k=1;
j=1;
for alpha=0:0.01:0.9;

    Ss=S;
    for i=2:N
        Ss(i)=Ss(i-1)+alpha*(S(i)-Ss(i-1));
    end

    Ssb=Ss;

    for i=N-1:-1:1
       Ssb(i)=Ssb(i+1)+alpha*(Ss(i)-Ssb(i+1));
    end

    %plot(T,Ssb,'g')

    [ Id_exp,Iv_exp ] = IdIv( S,Ssb );
    Id_expM(j)=Id_exp;
    Iv_expM(j)=Iv_exp;
    j=j+1;
    if Id_exp<Id_R && Iv_exp<Iv_R
        alphaM(k)=alpha;

        subplot(4,2,k)
        plot(T,Ssb,'r')
        hold on
        plot(T,S)
        title(alpha)
        k=k+1;
    end
end
