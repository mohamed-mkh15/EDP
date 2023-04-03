clear;
close all;
load('noisy_surface.mat')
load('true_surface.mat')

noisy=noisy_surface;
true=true_surface;

[x,y]=size(true);
N=x*y;

VarianceMatrix=(noisy-true).^2;

Variance=(1/(N-1))*sum(reshape(VarianceMatrix,[],1));

Z=noisy;
k=0;
for alpha=0.1:0.1:0.9;
    k=k+1;
    Zs1=Z;
    for j=1:y
        for i=2:x
            Zs1(i,j)=Zs1(i-1,j)+alpha*(Z(i,j)-Zs1(i-1,j));
        end
    end
    Zs2=Zs1;

    for j=1:y
        for i=x-1:-1:1
            Zs2(i,j)=Zs2(i+1,j)+alpha*(Zs1(i,j)-Zs2(i+1,j));
        end
    end

    Zs3=Zs2;

    for i=1:x
        for j=2:y
            Zs3(i,j)=Zs3(i,j-1)+alpha*(Zs2(i,j)-Zs3(i,j-1));
        end
    end

    Zs4=Zs3;

    for i=1:x
        for j=y-1:-1:1
            Zs4(i,j)=Zs4(i,j+1)+alpha*(Zs3(i,j)-Zs4(i,j+1));
        end
    end
    
    subplot(3,3,k)
    mesh(Zs4)
    title(alpha)
    VarianceMatrixS=(Zs4-true).^2;

    VarianceS=(1/(N-1))*sum(reshape(VarianceMatrixS,[],1));
    m(k) = VarianceS;
end
disp(m)