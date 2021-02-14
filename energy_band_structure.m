clear;
clc;


a0=2.46;
a1=a0*[0.5,sqrt(3)/2,0];
a2=a0*[-0.5,sqrt(3)/2,0];
a3=a0*[0,0,20];
N=1000;%--in unit of 1000
omega=dot(cross(a1,a2),a3);
b1=2*pi*cross(a2,a3)/omega;
b2=2*pi*cross(a3,a1)/omega;
b3=2*pi*cross(a1,a2)/omega;
E=eye(4.732*N,2);
G=[0,0];
K=[(b1(1)-b2(1))/3,(b1(2)-b2(2))/3];
M=[b1(1)/2,b1(2)/2];
%--k points from M to Gamma in reciprocal space
for k=1:+1:1.732*N
    E(k,1:2)=Hmatrix(M(1)+(G(1)-M(1))/(1.732*N)*k,M(2)+(G(2)-M(2))/(1.732*N)*k,1,2.8);
end
%--k points from Gamma to K in reciprocal space
for k=1.732*N+1:+1:3.732*N
    E(k,1:2)=Hmatrix(G(1)+(K(1)-G(1))/(2*N)*(k-1.732*N),G(2)+(K(2)-G(2))/(2*N)*(k-1.732*N),1,2.8);
end
%--k points from K to M in reciprocal space
for k=3.732*N+1:+1:4.732*N
    E(k,1:2)=Hmatrix(K(1)+(M(1)-K(1))/N*(k-3.732*N),K(2)+(M(2)-K(2))/N*(k-3.732*N),1,2.8);
end

figure (1)
for k=1:+1:4.732*N
    plot(k,E(k,1:2),'k.','MarkerSize',15);hold on;
end
%--show the boundary of the K points
xticks([0 1.732*N 3.732*N 4.732*N])
xticklabels({'M','г','K','M'})
line([1.732*N,1.732*N],[-8,12],'Color','k','LineWidth',1);hold on;
line([3.732*N,3.732*N],[-8,12],'Color','k','LineWidth',1);hold on;
axis([0 4.732*N -8 12]);
xlabel('K point distance');
ylabel('Eigenenergy(ev)');
title('E-k diagram with Tight binding (1NN) in graphene');

function Eeigen=Hmatrix(kx,ky,e,t) 
    a0=2.46;
    a1=a0*[0.5,sqrt(3)/2,0];
    a2=a0*[-0.5,sqrt(3)/2,0];
    %--d1-d3 is the vetor from B to A <A|H|B>
    d1=[-(a1(1)+a2(1))/3+a2(1),-(a1(2)+a2(2))/3+a2(2)];
    d2=[-(a1(1)+a2(1))/3+a1(1),-(a1(2)+a2(2))/3+a1(2)];
    d3=[-(a1(1)+a2(1))/3,-(a1(2)+a2(2))/3];
    %--d4-d6 is the vector from A to B <B|H|A>
    d4=-d3;
    d5=-d2;
    d6=-d1;
    H=zeros(2);
    H(1,1)=e;
    H(2,2)=e;
    H(1,2)=(exp(1i*dot([kx,ky],d1))+exp(1i*dot([kx,ky],d2))+exp(1i*dot([kx,ky],d3)))*t;
    H(2,1)=(exp(1i*dot([kx,ky],d4))+exp(1i*dot([kx,ky],d5))+exp(1i*dot([kx,ky],d6)))*t;
    Eeigen=eig(H);
end