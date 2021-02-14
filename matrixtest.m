clc
clear

a0=2.4669;
a1=a0*[0.5,sqrt(3)/2,0];
a2=a0*[-0.5,sqrt(3)/2,0];
a3=a0*[0,0,20.5];
omega=dot(cross(a1,a2),a3);
b1=2*pi*cross(a2,a3)/omega;
b2=2*pi*cross(a3,a1)/omega;
b3=2*pi*cross(a1,a2)/omega;
K=[(b1(1)-b2(1))/3,(b1(2)-b2(2))/3];
V=Hmatrix(K(1),K(2));
V
%construct the H(k)matrix with the 3-Nearest-Neighbor 
function [Eigenvector]=Hmatrix(kx,ky) 
    data=importdata("wannier90_hr_to_3NN.dat");
    dmatrix=ones(26,4);
    for m =1:length(data(:,1))
        [dmatrix(m,1),dmatrix(m,2),dmatrix(m,3),dmatrix(m,4)]=getdlengthandvector(data(m,1),data(m,2),data(m,3),data(m,4),data(m,5));
    end
    H=zeros(2,2);
    for m =1:+1:length(data(:,1))
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
    [V,~]=eig(H);
    Eigenvector=V;
end
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