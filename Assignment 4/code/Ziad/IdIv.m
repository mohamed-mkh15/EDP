function [ Id,Iv ] = IdIv( Z,X )
Id=sum((Z-X).^2);

xi=X(1:end-2);
for i=1:length(X)-2
  xi(i)=(X(i+2)-2*X(i+1)+X(i))^2;  
end
Iv=sum(xi);
end

