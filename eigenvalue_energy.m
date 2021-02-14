clear;
clc;

N=2;
output=zeros(N);
t=-0.6;
E=1;
for n=1:1:N
    H=zeros(n);
    for i=1:1:n
        for j=1:1:n
            if  j==i
                H(i,j)=E;
            elseif (i>=1)&&(j==i+1)||(j==i-1)
                H(i,j)=t;
            end
        end
    end
    [eigenvector,eigenvalue]=eig(H);
    for k=1:1:n
       output(k,n)=eigenvalue(k,k); 
    end
end
figure(1)
for i=1:1:N
   for n=1:1:i
      plot(i,output(n,i),'k.','MarkerSize',15);hold on; 
   end
end
figure(2)
for i=1:1:N
   A(i)=output(i,N);
   B=sort(A);
end
for i=1:1:N
    plot(i,B(i),'k.','MarkerSize',15);hold on; 
end
figure(3) %figure(1)energy/t
for i=1:1:N
   for n=1:1:i
      plot(i,(output(n,i)-E)/abs(t),'k.','MarkerSize',15);hold on; 
   end
end



