%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forest-fire spreading simulation using a default probability of burning
% for each tree and a default probability of staying on fire duing each
% iteration.
%
%
% Made by: Oscar A. Nieves
% Made in: 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;

% Inputs
n = 100; % trees per row and column
xvec = linspace(-5,5,n);
yvec = xvec;
forest = rand(n,n);
density = 0.8;
forest(forest > 1-density) = 1;
forest(forest <= 1-density) = 0;
vec = find(forest == 1);
index1 = randi(length(vec));

% Set fire position from available trees
forest(index1) = 0.5;

% Set square around forest (closed forest)
forest(1,:) = 0;
forest(end,:) = 0;
forest(:,1) = 0;
forest(:,end) = 0;

% Plot forest
figure(1);
imagesc(forest); % show forest at beginning

%% Percolation
iterations = 200;
kmax = sum(forest);
tree_prob = 0.9;
probability_off = 0.05;

% Apply mask (Select shape of forest and land)
mask1 = ones(n,n);
[XL,YL] = meshgrid(xvec,yvec);
R = ceil(n/2);
mask1(XL.^2 + YL.^2 <= 3.^2 & XL.^2 + YL.^2 >= 2.^2 ) = 0;
forest = forest.*mask1;
forest1 = forest;
plot_forest(:,:,1) = forest;

for ii = 1:iterations-1
    
    % Find trees on fire across entire domain
    [fx,fy] = find(forest == 0.5);
    
    % Apply probabilities of burning based on current trees on fire and
    % those adjacent to them
    if fx ~= 1
        for kk = 1:length(fx)
            if rand(1) < tree_prob && forest(fx(kk)-1,fy(kk)) ~= 0
                forest(fx(kk)-1,fy(kk)) = 0.5;
            end
        end
    end
    if fx ~= n
        for kk = 1:length(fx)
            if rand(1) < tree_prob && forest(fx(kk)+1,fy(kk)) ~= 0
                forest(fx(kk)+1,fy(kk)) = 0.5;
            end
        end
    end
    if fy ~= 1
        for kk = 1:length(fy)
            if rand(1) < tree_prob && forest(fx(kk),fy(kk)-1) ~= 0
                forest(fx(kk),fy(kk)-1) = 0.5;
            end
        end
    end
    if fy ~= n
        for kk = 1:length(fy)
            if rand(1) < tree_prob && forest(fx(kk),fy(kk)+1) ~= 0
                forest(fx(kk),fy(kk)+1) = 0.5;
            end
        end
    end
    
    % Update value of the forest plot array
    plot_forest(:,:,ii+1) = forest;
    
    % Update burning trees
    R = rand(n,n);
    forest(forest1==0.5 & forest==0.5 & R<probability_off) = 0;
    forest1 = forest;
end

%% Plot
close all;
figure(2);
set(gcf,'color','w');
for ii = 1:iterations
    imagesc(plot_forest(:,:,ii));
    colormap([0 0 0; 1 0 0; 0 1 0]);
    xlabel('x'); ylabel('y');
    title('Fire percolation');
    set(gca,'FontSize',20);
    axis tight;
    drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%