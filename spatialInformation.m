%Written by Guillaume Etter and edited by Tori-Lynn 

function info_per_frame = spatialInformation(ms,behav)

min_speed = 2;
bin_size = 3;

% Binning euclidian space
X_bin_vector = min(behav.position(:,1)):bin_size:max(behav.position(:,1));
Y_bin_vector = min(behav.position(:,2)):bin_size:max(behav.position(:,2));

%% Extracting data from structure
calcium_time = ms.time;
dt = mode(diff(calcium_time));
Fs = 1/dt;
    %% Re-aligning calcium and behavior times
    [unique_behavtime, IAbehav, ICbehav]=unique(behav.time);
    [unique_mstime, IAms, ICms]=unique(ms.time);

    % Interpolate the position of the mouse on the linear track
    interpolated_X=interp1(behav.time(IAbehav),behav.position(IAbehav,1),ms.time(IAms));
    interpolated_Y=interp1(behav.time(IAbehav),behav.position(IAbehav,2),ms.time(IAms));
    interpolated_speed=interp1(behav.time(IAbehav),behav.speed(IAbehav),ms.time(IAms));
    % Calculating spatial firing
    nbin = length(X_bin_vector)*length(Y_bin_vector); % Used for MI calculation
    I_frame = zeros(1,length(ms.Binary(1,:)));
    MI = zeros(1,length(ms.Binary(1,:)));
    firing_matrix = zeros(length(Y_bin_vector),length(X_bin_vector));
    
for cell_ID = 1:length(ms.Binary(1,:))
    calcium_trace = ms.RawTraces(:,cell_ID);
    binarized_trace = ms.Binary(:,cell_ID);
    avg_firing_prob(1,cell_ID) = sum(binarized_trace)./length(ms.time); %Expressed in probability of firing (<1)
    current_bin = 0;
    linear_position_vector = 0*ms.time;
    binarized_spatial_vector = 0*ms.time;
    

    for j = 1:length(Y_bin_vector)
        for i = 1:length(X_bin_vector)
            current_bin = current_bin+1;
            binarized_spatial_vector = 0*ms.time;
            position_idx = find(interpolated_X>X_bin_vector(i)-bin_size & interpolated_X < X_bin_vector(i) & interpolated_Y>Y_bin_vector(j)-bin_size & interpolated_Y < Y_bin_vector(j));

            if ~isempty(position_idx)
                position_idx(interpolated_speed(position_idx)<min_speed)=[]; % Remove indices for particular speeds

                linear_position_vector(position_idx) = current_bin; %Later used for MI calculation

                binarized_spatial_vector(position_idx)=1;
                prob_position_bin = length(position_idx)./length(ms.time);

                spikes_in_bin = sum(binarized_trace(position_idx),'omitnan');
                firing_probability = spikes_in_bin./length(position_idx);

                firing_in_bin_idx = find(binarized_trace == 1 & binarized_spatial_vector == 1);
                silent_in_bin_idx = find(binarized_trace == 0 & binarized_spatial_vector == 1);

                joint_prob_firing_in_bin = length(firing_in_bin_idx)/length(ms.time);
                joint_prob_silent_in_bin = length(silent_in_bin_idx)/length(ms.time);

                I_sec_temp = prob_position_bin*firing_probability*log2(firing_probability+eps/avg_firing_prob(1,cell_ID));
                if isnan(I_sec_temp)
                    I_sec_temp = 0;
                end
                I_frame(1,cell_ID) = I_frame(1,cell_ID) + I_sec_temp;

                MI_temp = joint_prob_firing_in_bin*log2(joint_prob_firing_in_bin+eps/(prob_position_bin*firing_probability)) + joint_prob_silent_in_bin*log2(joint_prob_silent_in_bin+eps/(prob_position_bin*(1-firing_probability)));

                if isnan(MI_temp);
                    MI_temp = 0;
                end
                MI(1,cell_ID) = MI(1,cell_ID) + MI_temp;

                firing_matrix(cell_ID,j,i) = joint_prob_firing_in_bin;

            else
                firing_matrix(cell_ID,j,i) = 0;
            end
        end
    end

end 

cellNum = 3; 


    color_vector = ms.time/ms.time(end);
    color_vector(:,3) = color_vector(:,1);
    color_vector(:,2) = flipud(color_vector(:,1));
    
    % Converting back to pixels
    interpolated_X_pix = interpolated_X*behav.ROI(3)/behav.trackLength;
    interpolated_Y_pix = interpolated_Y*behav.ROI(3)/behav.trackLength;
    
    %% Spatial firing
    clf
    plotting_fig = gcf;
    set(plotting_fig,'Name',strcat('Cell ID: ',num2str(cellNum)),'NumberTitle','off')
    subplot(3,2,[1 3]);
    imshow(3*behav.background); hold on
    plot(interpolated_X_pix,interpolated_Y_pix,'color',[1 1 1 0.3]); hold on;
    scatter(interpolated_X_pix(logical(ms.Binary(:,cellNum))),interpolated_Y_pix(logical(ms.Binary(:,cellNum))), [100],color_vector(logical(ms.Binary(:,cellNum)),:), '.');
    ax1=gca;
    colormap(ax1,color_vector);
    %cb=colorbar(ax1,'southoutside','Ticks',[0,1],'TickLabels',{'Start','End'});
    %cb.Label.String = 'Time';
    daspect([1 1 1])
    ax1.YDir = 'Reverse';
    set(ax1, 'Visible','off');hold off
    
    subplot(3,2,[2 4]);
    spatial_firing_heatmap = imagesc(squeeze(firing_matrix(cellNum,:,:)));
    ax=gca;
    ax.CLim = [0 0.01];
    hcb=colorbar;
    title(hcb,'Firing/occupancy joint probability')
    colormap 'jet';
    shading interp;
    set (gca,'DataAspectRatio',[1 1 1],'YDir','Reverse',...
        'XTick',[],...
        'YTick',[]);
    hold on;
    mouse_path = plot((interpolated_X+bin_size)./bin_size,(interpolated_Y+bin_size)./bin_size, 'color', [1 1 1 0.3]);
    hold off
    
    subplot(3,2,[5 6]);
    plot(ms.time/1000,ms.RawTraces(:,cellNum),'color','black');
    hold on
    plot(ms.time/1000,ms.Binary(:,cellNum),'color','red');
    ax7=gca;
    ax7.XLim= [0 ms.time(end)/1000];
    ax7.YLim= [0 max(ms.RawTraces(:,cellNum))];
    colormap(ax7,color_vector);
    cb=colorbar(ax7,'southoutside','Ticks',[0,1],'TickLabels',{'Start','End'});
    
    title(strcat('Bits/Frame: ',num2str(I_frame(1,cellNum)),', MI: ', num2str(MI(1,cellNum))));
    drawnow


info_per_frame = I_frame; 

end 