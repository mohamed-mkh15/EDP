clc 
clear 
n=500;
m=1; 
% Crating Arrays to save the values of Errors over 500 runs. 
XF_Err=ones(500,n); XF_E=ones(1,n);
YF_Err=ones(500,n); YF_E=ones(1,n);
DF_Err=ones(500,n); DF_E=ones(1,n);
BF_Err=ones(500,n); BF_E=ones(1,n); 
XP_Err=ones(500,n); XP_E=ones(1,n);
YP_Err=ones(500,n); YP_E=ones(1,n);
XE_Err=ones(500,n-m); XE_E=ones(1,n-m);
YE_Err=ones(500,n-m); YE_E=ones(1,n-m);
DE_Err=ones(500,n-m); DE_E=ones(1,n-m);
BE_Err=ones(500,n-m); BE_E=ones(1,n-m);
for j=1:500
    N=500;
    T=2;
    sigma_D=50;
    sigma_B=0.004;
    sigma_Ba=0.001;
    sigma_a=0.3;
    x=ones(1,N); x(1)=1000;
    y=ones(1,N); y(1)=1000;
    vx=ones(1,N); vx(1)=100;
    vy=ones(1,N); vy(1)=100;
    ax=normrnd(0,sigma_a,1,N);
    ay=normrnd(0,sigma_a,1,N);
    D=ones(1,N); D(1)=sqrt(x(1)^2+y(1)^2);
    B=ones(1,N); B(1)=atan2(x(1),y(1));
    Dm=ones(1,N); Dm(1)=D(1)+normrnd(0,sigma_D);
    Bm=ones(1,N); Bm(1)=B(1)+normrnd(0,sigma_B);
    for i=2:N
        x(i)=x(i-1)+vx(i-1)*T+(ax(i-1)*T^2)/2;
        y(i)=y(i-1)+vy(i-1)*T+(ay(i-1)*T^2)/2;
        vx(i)=vx(i-1)+ax(i-1)*T;
        vy(i)=vy(i-1)+ay(i-1)*T;
        D(i)=sqrt(x(i)^2+y(i)^2);
        B(i)=atan2(x(i),y(i));
        if floor(i/2)==i/2
            Bm(i)=B(i)+normrnd(0,sigma_Ba);
            Dm(i)=NaN;
        else
            Dm(i)=D(i)+normrnd(0,sigma_D);
            Bm(i)=B(i)+normrnd(0,sigma_B);
        end 
    end
    XF=ones(4,N);
    XE=ones(4,N-1);
    DE=ones(1,N-1); DF=ones(1,N);
    BE=ones(1,N-1); BF=ones(1,N);
    X0=[Dm(1)*sin(Bm(1));
       (Dm(3)*sin(Bm(3))-Dm(1)*sin(Bm(1)))/(2*T);
        Dm(1)*cos(Bm(1));
       (Dm(3)*cos(Bm(3))-Dm(1)*cos(Bm(1)))/(2*T)];
    P=diag([10^4 10^4 10^4 10^4 ]);
    phi=[1 T 0 0;
        0 1 0 0;
        0 0 1 T;
        0 0 0 1];
    G=[(T^2)/2 0;
        T 0;
        0 (T^2)/2;
        0 T];
    Q=G*G'.*sigma_a^2;
    dHa=getdH2(X0);
    Ra=sigma_Ba^2;
    Ka=P*dHa'/(dHa*P*dHa'+Ra);
    R=[sigma_D^2 0;
        0 sigma_B^2];
    dH=getdH(X0);
    K=P*dH'/(dH*P*dH'+R);
    for i=4:N
        if floor(i/2)==i/2
            Z=Bm(i);
            Xp=phi*X0;
            Ha=getHH(Xp);
            Xf=Xp+Ka*(Z-Ha);
            XF(:,i)=Xf;
            DF(i)=sqrt(Xf(1)^2+Xf(3)^2);
            BF(i)=atan2(Xf(1),Xf(3));
            Xe=phi*Xf;
            XE(:,i)=Xe;
            DE(i)=sqrt(Xe(1)^2+Xe(3)^2);
            BE(i)=atan2(Xe(1),Xe(3));
            Pex=phi*P*phi'+Q;
            dHa=getdH2(Xp);
            Ka=Pex*dHa'/(dHa*Pex*dHa'+Ra);
            P=(eye(4)-Ka*dHa)*Pex;
            X0=Xf;
        else
            Z=[Dm(i);
                Bm(i)];
            Xp=phi*X0;
            H=getH(Xp);
            Xf=Xp+K*(Z-H);
            XF(:,i)=Xf;
            DF(i)=sqrt(Xf(1)^2+Xf(3)^2);
            BF(i)=atan2(Xf(1),Xf(3));
            Xe=phi*Xf;
            XE(:,i)=Xe;
            DE(i)=sqrt(Xe(1)^2+Xe(3)^2);
            BE(i)=atan2(Xe(1),Xe(3));
            Pex=phi*P*phi'+Q;
            dH=getdH(Xp);
            K=Pex*dH'/(dH*Pex*dH'+R);
            P=(eye(4)-K*dH)*Pex;
            X0=Xf;    
        end
        % Errors Calculation. 
        XF_Err(j,i)=(x(i)-Xf(1))^2;
        YF_Err(j,i)=(y(i)-Xf(3))^2; 
        XP_Err(j,i)=sqrt(P(1,1)); 
        YP_Err(j,i)=sqrt(P(3,3));
        DF_Err(j,i)=(D(i)-DF(i))^2;
        BF_Err(j,i)=(B(i)-BF(i))^2;
        if i<N
            XE_Err(j,i)=(x(i+1)-Xe(1))^2;
            YE_Err(j,i)=(y(i+1)-Xe(3))^2;
            DE_Err(j,i)=(D(i+1)-DE(i))^2;
            BE_Err(j,i)=(B(i+1)-BE(i))^2;
        end
            
        
    end 

end 
for i=4:N
   XF_E(i)=sqrt((1/499)*sum(XF_Err(:,i))); 
   YF_E(i)=sqrt((1/499)*sum(YF_Err(:,i))); 
   DF_E(i)=sqrt((1/499)*sum(DF_Err(:,i))); 
   BF_E(i)=sqrt((1/499)*sum(BF_Err(:,i))); 
   XP_E(i)=(1/499)*sum(XP_Err(:,i)); 
   YP_E(i)=(1/499)*sum(YP_Err(:,i)); 
   if i<N
       XE_E(i)=sqrt((1/499)*sum(XE_Err(:,i))); 
       YE_E(i)=sqrt((1/499)*sum(YE_Err(:,i))); 
       DE_E(i)=sqrt((1/499)*sum(DE_Err(:,i))); 
       BE_E(i)=sqrt((1/499)*sum(BE_Err(:,i))); 
   end
end     

figure()
%Trajectory plotting for visualization
polarplot(B,D)
hold on
I = ~isnan(Bm) & ~isnan(Dm);
polarplot(Bm(I),Dm(I),'r')
polarplot(BF,DF,'g')
polarplot(BE,DE,'m')

title('Trajectory polt using Range and azimuth')
legend({'True trajectory','Measurements (Observer 1)','Filteration estimated','Extrapolation'},...
    'location','southwest')

%plotting errors
figure()
subplot(1,2,1)
plot(5:N,DF_E(5:N),'r'); hold on
plot(5:N-1,DE_E(5:N-1),'b'); hold on
plot(1:N,sigma_D*ones(1,N),'k')
legend('filteration True Estimation error','Extrapolation True estimation error','Standard deviation')
set(gcf,'position',[0,0,900,800]);
xlabel('points')
ylabel('value')
grid on

subplot(1,2,2)
plot(5:N,BF_E(5:N),'r'); hold on
plot(5:N-1,BE_E(5:N-1),'b'); hold on
plot(1:N,sigma_B*ones(1,N),'k'); hold on
plot(sigma_Ba*ones(1,N),'g')
title('True Estimation Error vs. Extrapolation error vs. standard deviation')
xlabel('points')
ylabel('value')
grid on

function [H]=getH(X)
 H=[sqrt(X(1)^2+X(3)^2); 
    atan2(X(1),X(3))];
end
function H=getHH(X)
 H=atan2(X(1),X(3));
end 
 function [dH]=getdH(X)
 dH=[X(1)/sqrt(X(1)^2+X(3)^2) 0 X(3)/sqrt(X(1)^2+X(3)^2) 0;
     X(3)/(X(1)^2+X(3)^2) 0 -X(1)/(X(1)^2+X(3)^2) 0];
 end
 function [dH]=getdH2(X)
 dH=[X(3)/(X(1)^2+X(3)^2) 0 -X(1)/(X(1)^2+X(3)^2) 0];
 end
 
 
