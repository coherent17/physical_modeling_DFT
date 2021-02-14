clc;
clear;
close all;

data1=importdata("k=( 1.698, 0.000, 0.000)_1BS_n1.dat");
data2=importdata("k=( 1.698, 0.000, 0.000)_2BS_n1.dat");
data3=importdata("k=( 1.698, 0.000, 0.000)_1BS_n2.dat");
data4=importdata("k=( 1.698, 0.000, 0.000)_2BS_n2.dat");
z = 10.591667;
order1 = find(data1.data(:,3)==z);

new_data(:,1) =  data1.data(order1,1);
new_data(:,2) =  data1.data(order1,2);
new_data(:,3) =  data1.data(order1,3);
new_data(:,4) =  data1.data(order1,4);
new_data(:,5) =  data1.data(order1,5);
% new_data(:,6) =  new_data(:,4).^2+new_data(:,5).^2;
new_data(:,6) =  new_data(:,4);


figure(1)
scatter(new_data(1:length(order1),1),new_data(1:length(order1),2),5,new_data(1:length(order1),6));hold on;
colorbar
title(colorbar,'...');
plot_carbon_atoms;
title(['Bloch sum (n=1) at K point on the surface z=' num2str(z) '(Å)']);
xlabel('x(Å)');
ylabel('y(Å)');
axis([-8 8 -8 8]);

function plot_carbon_atoms
    a0=2.4669;
    a1=a0*[0.5,sqrt(3)/2,0];
    a2=a0*[-0.5,sqrt(3)/2,0];
    X1(1)=(1/3)*(a1(1)+a2(1));
    Y1(1)=(1/3)*(a1(2)+a2(2));
    X2(1)=(2/3)*(a1(1)+a2(1));
    Y2(1)=(2/3)*(a1(2)+a2(2));
    X3(1)=a1(1);
    Y3(1)=a1(2);
    X4(1)=X3(1)+a2(1);
    Y4(1)=Y3(1)+a1(2);
    for i=1:+1:5
        theta=(60*i);
        X1(i+1)=X1(1)*cosd(theta)-Y1(1)*sind(theta);
        Y1(i+1)=Y1(1)*cosd(theta)+X1(1)*sind(theta);
        X2(i+1)=X2(1)*cosd(theta)-Y2(1)*sind(theta);
        Y2(i+1)=Y2(1)*cosd(theta)+X2(1)*sind(theta);
        X3(i+1)=X3(1)*cosd(theta)-Y3(1)*sind(theta);
        Y3(i+1)=Y3(1)*cosd(theta)+X3(1)*sind(theta);
        X4(i+1)=X4(1)*cosd(theta)-Y4(1)*sind(theta);
        Y4(i+1)=Y4(1)*cosd(theta)+X4(1)*sind(theta);
    end
    plot(X1,Y1,'k.','MarkerSize',20);hold on;
    plot(X2,Y2,'k.','MarkerSize',20);hold on;
    plot(0,0,'r.','MarkerSize',10);hold on;
    plot(X3,Y3,'r.','MarkerSize',10);hold on;
    plot(X4,Y4,'r.','MarkerSize',10);hold on;
end

