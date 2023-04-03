clear;
close all;
load('noisy_surface.mat')
load('true_surface.mat')

noisy=noisy_surface;
true=true_surface;

figure

mesh(true)

[x,y]=size(true);
N=x*y;

VarianceMatrix=(noisy-true).^2;

Variance=(1/(N-1))*sum(reshape(VarianceMatrix,[],1));

Z=noisy;
Zs1=Z;
alpha=0.335;
for j=1:y
    for i=2:x
        Zs1(i,j)=Zs1(i-1,j)+alpha*(Z(i,j)-Zs1(i-1,j));
    end
end
Zs3=Zs1;

for i=1:x
    for j=2:y
        Zs3(i,j)=Zs3(i,j-1)+alpha*(Zs1(i,j)-Zs3(i,j-1));
    end
end

hold on
mesh(Zs3,5*Zs3)


VarianceMatrixS=(Zs3-true).^2;

VarianceS=(1/(N-1))*sum(reshape(VarianceMatrixS,[],1));