function [Max_Value,Min_Value,T,X] = N_density(begin_year,xx,yy,zz,X)

[N,T] = size(X);

Year_0 = begin_year; % begin_year
Year_1 = Year_0+T-1;

%% Calculate the kernel density estimates for each year and plot the graphs.
figure(1)

for i = 1:T
     f(:,i) = ksdensity(X(:,i));
     plot(f)
     hold on
end

j = 0;
for i = Year_0:Year_1
      j = j+1;
      str{j} = char([num2str(i),'年']);
end

legend(str)

f = f';

%% Determine the X-axis and Y-axis ticks for the 3D kernel density estimate plot.
Min_Value = min(min(X));
Max_Value = max(max(X));

M = length(f);

x = linspace(Min_Value,Max_Value,M);

y = Year_0:Year_1;

[x,y] = meshgrid(x,y); 

%% Three-dimensional kernel density estimate result plot.
figure(2) 
mesh(x,y,f)
xlabel(xx);
ylabel(yy);
zlabel(zz);

set(gca,'FontName','Times New Roman','FontSize',30);
xlabel('\fontname{Times New Roman}\fontsize{37}CGPE');
ylabel('\fontname{Times New Roman}\fontsize{37}Year');
zlabel('\fontname{Times New Roman}\fontsize{37}Kernel estimates');
%title('\fontname{Times New Roman}\fontsize{44}(a) The whole study area')
%title('\fontname{Times New Roman}\fontsize{44}(b) Shanxi')
%title('\fontname{Times New Roman}\fontsize{44}(c) Henan')
%title('\fontname{Times New Roman}\fontsize{44}(d) Inner Mongolia')
%title('\fontname{Times New Roman}\fontsize{44}(e) Shaanxi')
%title('\fontname{Times New Roman}\fontsize{44}(f) Gansu')
%title('\fontname{Times New Roman}\fontsize{44}(g) Qinghai')
%title('\fontname{Times New Roman}\fontsize{44}(h) Ningxia')
end

