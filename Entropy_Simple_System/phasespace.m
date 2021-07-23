function [ Y ] = phasespace(signal,dim,tau)


N = length(signal);
% Total points on phase space 
T=N-(dim-1)*tau;
% Initialize the phase space
Y=zeros(T,dim);

for i=1:T
   Y(i,:)= signal(i+(dim-1)*tau-sort((0:dim-1),'descend')*tau)';
end

sizeY=size(Y,2);

if nargout == 0
    if sizeY == 2
        plot(Y(:,1),Y(:,2));
        xlabel('y1','FontSize',10,'FontWeight','bold');
        ylabel('y2','FontSize',10,'FontWeight','bold');
        get(gcf,'CurrentAxes');
        set(gca,'FontSize',10,'FontWeight','bold');
        grid on;
    else
        plot3(Y(:,1),Y(:,2),Y(:,3));
        xlabel('y1','FontSize',10,'FontWeight','bold');
        ylabel('y2','FontSize',10,'FontWeight','bold');
        zlabel('y3','FontSize',10,'FontWeight','bold');
        get(gcf,'CurrentAxes');
        set(gca,'FontSize',10,'FontWeight','bold');
        grid on;
    end
end