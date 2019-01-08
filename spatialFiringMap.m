%{
    Here we create a spatial firing map showing the areas of the maze 
    in which the cells fired the most. 

Here we do not exclude cells that fire when the mouse is running less than
the 2cm/s threshold. Spatial information does however exclude these. 
%}

function spatialFiringMap(ms,behav)
if ~isfield(ms,'Binary')
    disp('Extracting binary information');
    ms = msExtractBinary(ms);
end
%Get the behavorial position of the mouse on the track:
[~, IAbehav, ~]=unique(behav.time);
[~, IAms, ~]=unique(ms.time);

interpolated_X=interp1(behav.time(IAbehav),behav.position(IAbehav,1),ms.time(IAms),'pchip');
interpolated_Y=interp1(behav.time(IAbehav),behav.position(IAbehav,2),ms.time(IAms),'pchip');

bin_sizes = 1; 
length_of_field_x = round(max(interpolated_X)/bin_sizes)*bin_sizes;
length_of_field_y = round(max(interpolated_Y)/bin_sizes)*bin_sizes;

bins_x = length_of_field_x/bin_sizes;  
bins_y = length_of_field_y/bin_sizes; 
cells = length(ms.Binary(1,:));

cell_bins = zeros(cells, bins_y,bins_x);
cell_bins_all = zeros(bins_y,bins_x); 
x_pos = interpolated_X;
y_pos = interpolated_Y;
binSizes = bin_sizes; 

for cellNum = 1: cells
    for i = 1 : bins_x   
        %goes through each of the bins in the y direction
        for j =1:bins_y
           frames = 0;
           averageFiringInBins = 0;    
           frameVectorX = find(x_pos>(i-1)*binSizes & x_pos<=i*binSizes);
           frameVectorY = find(y_pos>(j-1)*binSizes  & y_pos<=j*binSizes);    
           frames = intersect(frameVectorX, frameVectorY);     
           averageFiringInBins =  sum(sum(ms.Binary(frames,cellNum),2)/length(ms.FiltTraces(1,cellNum)))/length(frames);
          averageFiringInBins_allCells =  sum(sum(ms.Binary(frames,:),2)/length(ms.FiltTraces(1,:)))/length(frames);
           if isnan(averageFiringInBins)
               averageFiringInBins = 0; 
           end        
           cell_bins(cellNum,j,i) = averageFiringInBins;     
           cell_bins_all(j,i) = averageFiringInBins_allCells;
        end 
    end
end 

%Trajectory of mouse running around with cell firing: 
figure  
plot(interpolated_X,interpolated_Y)
set(gca,'Ydir', 'reverse')
title('Mouse Trajectory')
xlabel('Length of the field (cm)')
ylabel('Width of the field (cm)')
hold on 

%Firing of a specified cell over the trajectory. 
%"frames" finds the indices of the non binary values within binariziedTraces 
frame = find(ms.Binary(:,15)); 
RGB = [rand(1) rand(1) rand(1)]; 
scatter(interpolated_X(frame),interpolated_Y(frame),[],RGB);   
hold off


imagesc(cell_bins_all); 
shading interp;
set(gca,'TickDir','out')
colormap('jet');
colorbar; 
title('Heat Maps of Cells Firing')
xlabel('Length of the field (cm)')
ylabel('Width of the field (cm)')

end 