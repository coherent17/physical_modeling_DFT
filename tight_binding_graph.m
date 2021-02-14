clc;
clear;

figure(1)
plotcell(0,0);
plotcell(0,1);
plotcell(1,0);
plotcell(1,1);
plotcell(-1,0);
plotcell(0,-1);
plotcell(-1,-1);
plotcell(1,-1);
plotcell(-1,1);

axis([-8 8 -8 8])
axis off;

function plotcell(n1,n2)
    a0=2.4669;
    a1=a0*[0.5,sqrt(3)/2,0];
    a2=a0*[-0.5,sqrt(3)/2,0];
    %plot the middle point in the cell
    plot(n1*a1(1)+n2*a2(1),n1*a1(2)+n2*a2(2),'k.','Markersize',1);hold on;
    
%     %the lattice point of the cell
%     %right:
%     plot(n1*a1(1)+n2*a2(1)+a1(1),n1*a1(2)+n2*a2(2),'g.','Markersize',5);hold on;
%     %left:
%     plot(n1*a1(1)+n2*a2(1)+a2(1),n1*a1(2)+n2*a2(2),'g.','Markersize',5);hold on;
%     %up:
%     plot(n1*a1(1)+n2*a2(1),n1*a1(2)+n2*a2(2)+a1(2),'g.','Markersize',5);hold on;
%     %down:
%     plot(n1*a1(1)+n2*a2(1),n1*a1(2)+n2*a2(2)-a1(2),'g.','Markersize',5);hold on;
    
    %Drawline of the boundary of the cell
    %up to right
    DrawLine(n1*a1(1)+n2*a2(1),n1*a1(2)+n2*a2(2)+a1(2),n1*a1(1)+n2*a2(1)+a1(1),n1*a1(2)+n2*a2(2));
    %right to down
    DrawLine(n1*a1(1)+n2*a2(1)+a1(1),n1*a1(2)+n2*a2(2),n1*a1(1)+n2*a2(1),n1*a1(2)+n2*a2(2)-a1(2));
    %down to left
    DrawLine(n1*a1(1)+n2*a2(1),n1*a1(2)+n2*a2(2)-a1(2),n1*a1(1)+n2*a2(1)+a2(1),n1*a1(2)+n2*a2(2));
    %left to up
    DrawLine(n1*a1(1)+n2*a2(1)+a2(1),n1*a1(2)+n2*a2(2),n1*a1(1)+n2*a2(1),n1*a1(2)+n2*a2(2)+a1(2));
    
    
    %plot the carbon atoms in the cell
    plot(n1*a1(1)+n2*a2(1)+(1/6)*(a1(1)+a2(1)),n1*a1(2)+n2*a2(2)+(1/6)*(a1(2)+a2(2)),'r.','Markersize',15);hold on;%the upper carbon atom in the cell
    plot(n1*a1(1)+n2*a2(1)-(1/6)*(a1(1)+a2(1)),n1*a1(2)+n2*a2(2)-(1/6)*(a1(2)+a2(2)),'b.','Markersize',15);hold on;%the down carbon atom in the cell
end

function DrawLine(x1,y1,x2,y2)
    line([x1,x2],[y1,y2],'Color','r','LineWidth',1);hold on;
end