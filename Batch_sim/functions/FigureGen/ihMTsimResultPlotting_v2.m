function ihMTsimResultPlotting_v2 (x, y, z, xLabel, yLabel, zLabel, zLow, zHigh,typePlot, countourVal )

% V2 implements a few changes for the surface plots:
% - remove grid lines
% - add in grey points (via scatter3) to show the grid of parameters that were simulated
%   but didnt make the time restriction

% Use turbo instead of JET colour to avoid false edges


% typePlot -> choose between scatter, surf or line plots
if strcmp(typePlot, 'scatter')
    subplot(2,1,1);
    scatter3( x, y,  z, 40, z, 'filled')
        ylabel(xLabel, 'FontSize', 20, 'FontWeight', 'bold')
        xlabel(yLabel, 'FontSize', 20, 'FontWeight', 'bold')
        zlabel(zLabel, 'FontSize', 20, 'FontWeight', 'bold')
        ax = gca;
        ax.FontSize = 20; 
        zlim([zLow zHigh])
        set(gcf,'position',[10,100,600,800])
        view(80,15)

        subplot(2,1,2); 
        scatter3( x, y,  z, 40, z, 'filled')
        ylabel(xLabel, 'FontSize', 20, 'FontWeight', 'bold')
        xlabel(yLabel, 'FontSize', 20, 'FontWeight', 'bold')
        zlabel(zLabel, 'FontSize', 20, 'FontWeight', 'bold')
        ax = gca;
        ax.FontSize = 20; 
        zlim([zLow zHigh])
        view(0,90)
        
       %% 
elseif strcmp(typePlot, 'surf')
    xu = unique(x,'sorted');
    yu = unique(y,'sorted');
    
    [xx, yy] = meshgrid (xu,yu);
    zz = zeros(size (xx));
    
    % First convert 1D variables into meshgrid
    for i = 1:length(x)
           
            %find indicie for x and y values to insert into Z
            tmpx = x(i);
            tmpy = y(i);
            
            idX = xx(1,:) == tmpx;
            idY = yy(:,1) == tmpy;
            
            zz (idY, idX) = z(i);
    end
    
    %% Could be an easy way to get grey values, take the NaN values below and extract them here!
    greyZ = zeros(size(zz));
    greyZ( zz == 0) = 1;
    greyX = xx;
    greyY = yy;

    % Convert to 1D for scatter3
    greyX= greyX(:);
    greyY= greyY(:);
    greyZ= greyZ(:);
    % remove spots that we want to plot sim results for
    greyX(greyZ == 0) = [];
    greyY(greyZ == 0) = [];
    greyZ(greyZ == 0) = [];

       % For surf to work, set all other values to NaN
    zz(zz ==0) = NaN; 
    
    % JAN 7TH 2021, temp remove the first plot.
%          subplot(2,1,1);
%          surf( xx, yy,  zz, 'EdgeColor', 'interp', 'FaceColor','interp' , 'FaceAlpha',0.3)
%             colormap turbo
%             ylabel(xLabel, 'FontSize', 20, 'FontWeight', 'bold')
%             xlabel(yLabel, 'FontSize', 20, 'FontWeight', 'bold')
%             zlabel(zLabel, 'FontSize', 20, 'FontWeight', 'bold')
%             ax = gca;
%             ax.FontSize = 20; 
%             zlim([zLow zHigh])
%             set(gcf,'position',[10,100,600,800])
%             hold on 
%             scatter3( x, y,  z, 40, z, 'filled')
%             view(80,15)

            %subplot(2,1,2); 
            figure;
            scatter3(greyX, greyY,greyZ, 40, [0.7 0.7 0.7], 'filled')
            hold on

            surf( xx, yy,  zz,  'EdgeColor', 'interp', 'FaceColor','interp' , 'FaceAlpha',0.3,'EdgeColor','none')
            if countourVal
                contour(xx, yy,  zz, 0:14)
            end
            colormap turbo
            ylabel(xLabel, 'FontSize', 20, 'FontWeight', 'bold')
            xlabel(yLabel, 'FontSize', 20, 'FontWeight', 'bold')
            zlabel(zLabel, 'FontSize', 20, 'FontWeight', 'bold')
            ax = gca;
            ax.FontSize = 20; 
            zlim([zLow zHigh])
            caxis([zLow zHigh])
            hold on 
            scatter3( x, y,  z, 40, z, 'filled')
            view(0,90)
            colorbar
            grid off
            hold off

            %%     
elseif strcmp(typePlot, 'line')
        
        % First sort based X
        [x,sortIdx] = sort(x,'ascend');
        % sort B using the sorting index
        y = y(sortIdx);
        z = z(sortIdx);

        subplot(2,1,1);
        plot3( x, y,  z,'LineWidth',3)
            ylabel(xLabel, 'FontSize', 20, 'FontWeight', 'bold')
            xlabel(yLabel, 'FontSize', 20, 'FontWeight', 'bold')
            zlabel(zLabel, 'FontSize', 20, 'FontWeight', 'bold')
            ax = gca;
            grid on
            ax.FontSize = 20; 
            zlim([zLow zHigh])
            set(gcf,'position',[10,100,600,800])
            view(80,15)

            subplot(2,1,2); 
            plot3( x, y,  z,'LineWidth',3)
            ylabel(xLabel, 'FontSize', 20, 'FontWeight', 'bold')
            xlabel(yLabel, 'FontSize', 20, 'FontWeight', 'bold')
            zlabel(zLabel, 'FontSize', 20, 'FontWeight', 'bold')
            ax = gca;
            grid on
            ax.FontSize = 20; 
            zlim([zLow zHigh])
            view(0,90)

    
    
else
    error('Select scatter, surf, or line for typePlot variable')
end
    

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    
    