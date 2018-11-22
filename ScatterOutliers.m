function stats = ScatterOutliers(a,b,j,fc,cm)
%
% ScatterOutliers(a,b,[j=1000,fc=true,cm='cool'])
%
% Scatter plot on top of a contour plot of the Mahalanobis lines and a robust regression line.
%
% If j > 0, the Mahalanobis distances are bootstrapped (i.e. Shepherd's pi correlation).
% 
% If fc is true the contour plot is filled, otherwise only the lines are drawn.
% 
% The final argument cm defines the colour map to be used. Default is 'cool'.
% If this begins with '-' the colour map is inverted (e.g. '-hot').
%
% Requires the Statistics Toolbox.
%

% Default parameters
if nargin < 3
    j = 1000;       % Number of bootstraps
    fc = true;      % Fill the contour plot
    cm = 'cool';    % Use colour map cool
elseif nargin < 4
    fc = true;
    cm = 'cool';
elseif nargin < 5
    cm = 'cool';
end

% Determine axes dimensions
hold off; scatter(a,b);
ra = xlim;
rb = ylim;
% Calculate Shepherd's pi correlation
[r p ol] = Shepherd(a,b,j);
% Bootstrapped Mahalanobis distances of all points in a grid
[x y] = meshgrid(ra(1):range(ra)/100:ra(2),rb(1):range(rb)/100:rb(2));
if j > 0
    gm = bsmahal([x(:) y(:)],[a b],j);
else
    gm = mahal([x(:) y(:)],[a b]);
end
gm = reshape(gm, size(x));
% Plot grid points as contours
if fc 
    contourf(x,y,gm,6:6:300);
else
    contour(x,y,gm,6:6:300);
end
hold on
% Plot the actual data 
scatter(a,b,100,'k','filled');
% Outliers are shown as open circles
scatter(a(ol),b(ol),100,'k','markerfacecolor','w');
% Linear fit on outlier-cleaned data
f = polyfit(a(~ol),b(~ol),1);
line(ra,ra*f(1)+f(2),'linewidth',3,'color','k');
% Invert colour map if desired
if cm(1) == '-'
    cm = flipud(colormap(cm(2:end)));
end
colormap(cm);
stats = [r p];
axis square
set(gca,'fontsize',15);
if j == 0
    % Without bootstrapping this is just Spearman's rho with outlier removal
    rhostr = 'rho';
else
    rhostr = 'pi';
end
p = round_decs(p,4);
if p == 0
    pstr = '<0.0001';
else
    pstr = num2str(p);
end
title([rhostr '=' num2str(round_decs(r,2)) ', p=' pstr]);