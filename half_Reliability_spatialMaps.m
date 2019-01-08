%{
     How we know the reliability of a spatial map

We know how much each cell fires in each bin in the maze 
Need to chop up the location a cell fire in half, you can
calculate the correlation by 

Correlation coefficiate: 

              r = (Zx1Zy1 +  Zx2Zy2 + Zx3Zy3)/(1-3) 


Split half reliablity, compare the correlation between location of the cell
firing in the first half vs. the second half. 

%}
 

function half_Reliability_spatialMaps(ms,behav, cell_ID, cell_bins_all)

min_speed = 2;
bin_size = 3;


[~, IAbehav, ~]=unique(behav.time);
[~, IAms, ~]=unique(ms.time);
% Binning euclidian space
X_bin_vector = min(behav.position(:,1)):bin_size:max(behav.position(:,1));
Y_bin_vector = min(behav.position(:,2)):bin_size:max(behav.position(:,2));

% Interpolate the position of the mouse on the linear track
interpolated_X=interp1(behav.time(IAbehav),behav.position(IAbehav,1),ms.time(IAms));
interpolated_Y=interp1(behav.time(IAbehav),behav.position(IAbehav,2),ms.time(IAms));
interpolated_speed=interp1(behav.time(IAbehav),behav.speed(IAbehav),ms.time(IAms));

%getting image1 and image2: 
halfFrames = length(ms.time(:,1))/2; %half the  frames of a session
if rem(halfFrames,1)~=0 
    halfFrames = halfFrames +0.5; 
end 
%We need to do this for the first half of the session and compare for the
%secondL 

%first half:
ms.FirstTime = ms.time(1:halfFrames,:); 
posXFirst = interpolated_X(1:halfFrames,1);
posYFirst = interpolated_Y(1:halfFrames,1);
SpeedFirst = interpolated_speed(1:halfFrames,1);
binarized_trace_first = ms.Binary(1:halfFrames,cell_ID);

for j = 1:length(Y_bin_vector);
    for i = 1:length(X_bin_vector);
        binarized_spatial_vector_first = 0*ms.FirstTime;
        position_idx_first = find(posXFirst>X_bin_vector(i)-bin_size & posXFirst < X_bin_vector(i) & posYFirst>Y_bin_vector(j)-bin_size & posYFirst < Y_bin_vector(j));
        
        if ~isempty(position_idx_first);
            position_idx_first(SpeedFirst(position_idx_first)<min_speed)=[]; % Remove indices for particular speeds            
            binarized_spatial_vector_first(position_idx_first)=1;
            firing_in_bin_idx_first = find(binarized_trace_first == 1 & binarized_spatial_vector_first == 1);            
            joint_prob_firing_in_bin_first = length(firing_in_bin_idx_first)/length(ms.FirstTime);            
            firing_matrix_first(j,i) = joint_prob_firing_in_bin_first;            
        else
            firing_matrix_first(j,i) = 0;
        end
    end
end

%Second half: 
ms.SecondTime = ms.time(halfFrames:end,:); 
posXSecond = interpolated_X(halfFrames:end,1);
posYSecond = interpolated_Y(halfFrames:end,1);
SpeedSecond = interpolated_speed(halfFrames:end,1);
binarized_trace_Second = ms.Binary(halfFrames:end,cell_ID);

for j = 1:length(Y_bin_vector);
    for i = 1:length(X_bin_vector);
        binarized_spatial_vector_Second = 0*ms.SecondTime;
        position_idx_Second = find(posXSecond>X_bin_vector(i)-bin_size & posXSecond < X_bin_vector(i) & posYSecond>Y_bin_vector(j)-bin_size & posYSecond < Y_bin_vector(j));
        
        if ~isempty(position_idx_Second);
            position_idx_Second(SpeedSecond(position_idx_Second)<min_speed)=[]; % Remove indices for particular speeds            
            binarized_spatial_vector_Second(position_idx_Second)=1;
            firing_in_bin_idx_Second = find(binarized_trace_Second == 1 & binarized_spatial_vector_Second == 1);            
            joint_prob_firing_in_bin_Second = length(firing_in_bin_idx_Second)/length(ms.SecondTime);            
            firing_matrix_Second(j,i) = joint_prob_firing_in_bin_Second;            
        else
            firing_matrix_Second(j,i) = 0;
        end
    end
end
  
%%----comparing the first half of the frames with the second: 
correlation_halfhalf = corr2(firing_matrix_first,firing_matrix_Second); 
if isnan(correlation_halfhalf)
    correlation_halfhalf =0; 
end 

%plotting a matrix that compares the first half with the second half: 

subplot(2,2,1);
imagesc(firing_matrix_first); 
shading interp;
set(gca,'TickDir','out')
colormap('jet');
colorbar; 
title('Heat Maps of Cells Firing in first half')
xlabel('Length of the field (cm)')
ylabel('Width of the field (cm)')

subplot(2,2,2);
imagesc(firing_matrix_Second); 
shading interp;
set(gca,'TickDir','out')
colormap('jet');
colorbar; 
title('Heat Maps of Cells Firing in second half')
xlabel('Length of the field (cm)')
ylabel('Width of the field (cm)')

subplot(2,2,[3 4]);
imagesc(cell_bins_all); 
shading interp;
set(gca,'TickDir','out')
colormap('jet');
colorbar; 
title(strcat('Correlation of half split reliability: ', num2str(correlation_halfhalf)))
xlabel('Length of the field (cm)')
ylabel('Width of the field (cm)')


end