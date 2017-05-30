

% Demo script for visualzing K-clusters (for K = 2)
% using Mulitdimensional scaling


curves_m = [curves;means(1,1);means(2,1)];
l = length(curves_m);
cost_mat = zeros(l,l);

for i=1:l
    for j=1:l
        cost_mat(j,i) = curves_m{j,1} - curves_m{i,1};
    end
end

% diagonal elements set to zero
for i=1:l
    cost_mat(i,i) = 0;
end


% visualizing the results with multidimensional scaling
K = cost_mat;
red = 2;
y = mdscale(K,red,'criterion','metricsstress');

% plotting result with the illustrations
y = y*10; % coordinate scaling for better visualization

figure;
for i=1:l
    test = curves_m{i,1}.points;
    test = test + repmat(y(i,1:2),length(test),1);
    if i== (l-1)
        line(test(:,1),test(:,2),'LineWidth',2,'color','red');
        hold on;
        line([test(end,1);test(1,1)],[test(end,2);test(1,2)],'LineWidth',2,'color','red');
        hold on;
    elseif i== l
        line(test(:,1),test(:,2),'LineWidth',2,'color','blue');
        hold on;
        line([test(end,1);test(1,1)],[test(end,2);test(1,2)],'LineWidth',2,'color','blue');
        hold on;
    else
        line(test(:,1),test(:,2),'LineWidth',2,'color','black');
        hold on;
        line([test(end,1);test(1,1)],[test(end,2);test(1,2)],'LineWidth',2,'color','black');
        hold on;
    end
    
    % member data indicator
    m =  mean(test);
    if any(i==elements{1,1})
        scatter(m(1,1),m(1,2),'filled','red'); hold on;
    elseif any(i==elements{2,1})
        scatter(m(1,1),m(1,2),'filled','blue'); hold on;
    end
end
axis square;
box on;
hold off;


