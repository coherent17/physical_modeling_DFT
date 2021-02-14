clc;
clear;

a0=2.4669/2;
a1=a0*[0.5,sqrt(3)/2,0];
a2=a0*[-0.5,sqrt(3)/2,0];

%the carbon atoms in the unit cell
X1(1)=(a1(1)+a2(1))/3;
Y1(1)=(a1(2)+a2(2))/3;
X1(2)=(a1(1)+a2(1))/3*2;
Y1(2)=(a1(2)+a2(2))/3*2;

X2(1)=a1(1)/3+4*a2(1)/3;
Y2(1)=a1(2)/3+4*a2(2)/3;
X2(2)=2*a1(1)/3+5*a2(1)/3;
Y2(2)=2*a1(2)/3+5*a2(2)/3;

X3(1)=4*a1(1)/3+a2(1)/3;
Y3(1)=4*a1(2)/3+a2(2)/3;
X3(2)=5*a1(1)/3+2*a2(1)/3;
Y3(2)=5*a1(2)/3+2*a2(2)/3;

X4(1)=4*a1(1)/3+4*a2(1)/3;
Y4(1)=4*a1(2)/3+4*a2(2)/3;
X4(2)=5*a1(1)/3+5*a2(1)/3;
Y4(2)=5*a1(2)/3+5*a2(2)/3;

%the lattice point to draw the boundary of the unit cell
L_X1=[0,a1(1),a1(1)+a2(1),a2(1)];
L_Y1=[0,a1(2),a1(2)+a2(2),a2(2)];

L_X2=L_X1+a2(1);
L_Y2=L_Y1+a2(2);

L_X3=L_X1+a1(1);
L_Y3=L_Y1+a1(2);

L_X4=L_X1+(a1(1)+a2(1));
L_Y4=L_Y1+(a1(2)+a2(2));

figure(1)
plot(X1,Y1,'k.','Markersize',10);hold on;
plot(X2,Y2,'k.','Markersize',10);hold on;
plot(X3,Y3,'k.','Markersize',10);hold on;
plot(X4,Y4,'k.','Markersize',10);hold on;
output_unit_cell(L_X1,L_Y1);
output_unit_cell(L_X2,L_Y2);
output_unit_cell(L_X3,L_Y3);
output_unit_cell(L_X4,L_Y4);
axis([-4.5 4.5 0 9]);
axis off;
title('Supercell graphene 2X2     unit(Å)');

function output_unit_cell(X,Y)
    plot(X,Y,'r.','Markersize',10);hold on;
    line([X(1),X(2)],[Y(1),Y(2)],'Color','r','LineWidth',1);hold on;
    line([X(2),X(3)],[Y(2),Y(3)],'Color','r','LineWidth',1);hold on;
    line([X(3),X(4)],[Y(3),Y(4)],'Color','r','LineWidth',1);hold on;
    line([X(4),X(1)],[Y(4),Y(1)],'Color','r','LineWidth',1);hold on;
end
