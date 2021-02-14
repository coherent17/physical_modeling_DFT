clc;
clear;
%lattice vector of the graphene in real space
a0=2.4669;
a1=a0*[0.5,sqrt(3)/2,0];
a2=a0*[-0.5,sqrt(3)/2,0];
a3=a0*[0,0,20.5];
%import the wannier 90 data with tight binding 3NN
data=importdata("wannier90_hr_to_3NN.dat");
%assign the matrix to store the d length and the d vector in x y z
dmatrix=ones(26-6,4);
for m =1:length(data(:,1))-6
    [dmatrix(m,1),dmatrix(m,2),dmatrix(m,3),dmatrix(m,4)]=getdlengthandvector(data(m,1),data(m,2),data(m,3),data(m,4),data(m,5));
end
N=1000;%--in unit of 1000
omega=dot(cross(a1,a2),a3);
b1=2*pi*cross(a2,a3)/omega;
b2=2*pi*cross(a3,a1)/omega;
b3=2*pi*cross(a1,a2)/omega;
E=ones(4.732*N,2);
G=[0,0];
K=[(b1(1)-b2(1))/3,(b1(2)-b2(2))/3];
M=[b1(1)/2,b1(2)/2];
%--k points from M to Gamma in reciprocal space
for k=1:+1:1.732*N
    E(k,1:2)=Hmatrix(M(1)+(G(1)-M(1))/(1.732*N)*k,M(2)+(G(2)-M(2))/(1.732*N)*k);
end
%--k points from Gamma to K in reciprocal space
for k=1.732*N+1:+1:3.732*N
    E(k,1:2)=Hmatrix(G(1)+(K(1)-G(1))/(2*N)*(k-1.732*N),G(2)+(K(2)-G(2))/(2*N)*(k-1.732*N));
end
%--k points from K to M in reciprocal space
for k=3.732*N+1:+1:4.732*N
    E(k,1:2)=Hmatrix(K(1)+(M(1)-K(1))/N*(k-3.732*N),K(2)+(M(2)-K(2))/N*(k-3.732*N));
end
figure (1)
for k=1:+1:4.732*N
    plot(k,real(E(k,1:2)+3.2),'m.','MarkerSize',15);hold on;
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

function [dlength,dx,dy,dz]=getdlengthandvector(nx,ny,nz,i,j)
    a0=2.4669;
    a1=a0*[0.5,sqrt(3)/2,0];
    a2=a0*[-0.5,sqrt(3)/2,0];
    a3=a0*[0,0,20.5];
    d_i=0*a1+0*a2+0*a3+(a1+a2)/3*i;
    d_j=nx*a1+ny*a2+nz*a3+(a1+a2)/3*j;
    d=d_i-d_j;
    dlength=sqrt(d(1)^2+d(2)^2+d(3)^2);
    dx=d(1);
    dy=d(2);
    dz=d(3);
end
%construct the H(k)matrix with the 3-Nearest-Neighbor 
function Eigenvalue=Hmatrix(kx,ky) 
    data=importdata("wannier90_hr_to_3NN.dat");
    dmatrix=ones(26-6-12,4);
    for m =1:length(data(:,1))-6-12
        [dmatrix(m,1),dmatrix(m,2),dmatrix(m,3),dmatrix(m,4)]=getdlengthandvector(data(m,1),data(m,2),data(m,3),data(m,4),data(m,5));
    end
    H=zeros(2,2);
    for m =1:+1:26-6-12
        i=data(m,4);
        j=data(m,5);
       if (i==1&&j==1)
           H(1,1)=H(1,1)+(exp(1i*dot([kx,ky],dmatrix(m,2:3))))*(data(m,6)+data(m,7)*1i);
       elseif (i==1&&j==2)
           H(1,2)=H(1,2)+(exp(1i*dot([kx,ky],dmatrix(m,2:3))))*(data(m,6)+data(m,7)*1i);
       elseif (i==2&&j==1)
           H(2,1)=H(2,1)+(exp(1i*dot([kx,ky],dmatrix(m,2:3))))*(data(m,6)+data(m,7)*1i);
       elseif (i==2&&j==2)
           H(2,2)=H(2,2)+(exp(1i*dot([kx,ky],dmatrix(m,2:3))))*(data(m,6)+data(m,7)*1i);
       end
    end
    Eigenvalue=eig(H);
end