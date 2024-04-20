
%% load task label (language task) 
flagSur = 0;
disp(['loading task label...'])    
if flagTask == 4 % Motor Task
    % load task label of each subject
    foldername = [main_folder,'/Sample Data/Motor Task/TaskLabel'];
    cd(foldername)
    name = dir(pwd) ;   
    file_name2 = ['MotorTaskLabelAllSubject.mat'];                    
    load (file_name2);     
    for isubject = 1:150   
        fullTime_allsubject{isubject} = TaskLabel_AllSubject_motor{isubject};  
    end
end


if flagTask == 4 % Motor Task
    % define parameters
    session_duration = 17;   %(currently considering this as the number of frames 12 seconds/0.72 )
    disp('[extracting task specific spiral centre distributions...]')    
%% MAIN SUBJECT FOR LOOP
    no_of_trial=0;
    for subject = 1:150
        disp('[Loading Raw Phase Signals...]')
        tic
        if hemisphere == 1
            % load spatiotemporally bandpass filtered (smoothed) real data
            folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
            cd(folder_name)
            filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_LEFT_sub',num2str(subject),'.mat'];
            load(filename)
                            
            DataIn_smooth = DataOut(1:175,:,1:283);
            % load temporally bandpass filtered (unsmoothed) real data
            folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
            cd(folder_name)
            filename = ['Preprocessed_temporalbandpass_data_motor_task_LEFT_sub',num2str(subject),'.mat'];
            load(filename)
            DataIn_unsmooth = DataOut(1:175,:,1:283);
        elseif hemisphere == 2
            % load spatiotemporally bandpass filtered (smoothed) real data
            folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
            cd(folder_name)
            filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_RIGHT_sub',num2str(subject),'.mat'];
            load(filename)
            
            DataIn_smooth = DataOut(1:175,:,1:283);
            % load temporally bandpass filtered (unsmoothed) real data
            folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
            cd(folder_name)
            filename = ['Preprocessed_temporalbandpass_data_motor_task_RIGHT_sub',num2str(subject),'.mat'];
            load(filename)
            DataIn_unsmooth = DataOut(1:175,:,1:283);
        end
        toc

        disp('[Loading Detected Spirals...]')
        tic
        foldername = [main_folder,'/Sample Data/Motor Task/Spiral Detected'];
        cd(foldername)
        if hemisphere == 1
            filename = ['Spiral_detected_surfilt_motor_task_LEFT_sub',num2str(subject),'.mat']; 

        elseif hemisphere == 2
            filename = ['Spiral_detected_surfilt_motor_task_RIGHT_sub',num2str(subject),'.mat']; 

        end
        load(filename)
        toc
        %  Find phase map of fMRI signal
        smooth_phase_map_real = nan(size(DataIn_smooth(1:175,:,:)));
        raw_phase_map_real = nan(size(DataIn_unsmooth(1:175,:,:)));
        for irow = 1:size(smooth_phase_map_real,1)
            for icol = 1:size(smooth_phase_map_real,2)
                temp1 = DataIn_smooth(irow,icol,:);
                if nansum(temp1(:))~=0
                smooth_phase_map_real(irow,icol,:) = angle(hilbert(temp1(:)));
                smooth_phase_map{subject} = smooth_phase_map_real;
                end
                temp1 = DataIn_unsmooth(irow,icol,:);
                if nansum(temp1(:))~=0
                raw_phase_map_real(irow,icol,:) = angle(hilbert(temp1(:)));  % frequency represented as rate of change in phase
                end        
            end
        end


        if flagSur == 0
            spiral_filt_nega = spiral_filt_nega_real_95perc_extend; % spiral full to calculate average pos and centre of mass
            spiral_filt_pos = spiral_filt_pos_real_95perc_extend;

            spiral_filt_nega_centre = spiral_filt_nega_real_centreONLY_95perc_extend; % Spiral centres to calculate distance
            spiral_filt_pos_centre = spiral_filt_pos_real_centreONLY_95perc_extend;
        end

        if flagRest == 0
            temp1_accu_1subject_spiral_pos = zeros(175,251,283);
            temp1_accu_1subject_spiral_neg = zeros(175,251,283);
        elseif flagRest == 1
            temp1_accu_1subject_spiral_pos = zeros(175,251,1199);   % Should see if rest data contains information on breaks between tasks or something
            temp1_accu_1subject_spiral_neg = zeros(175,251,1199);
        end

        % clockwise spirals
        for ipatt = 1:size(spiral_filt_nega,1)
            for t = 1:size(spiral_filt_nega,2)
                temp1 = full(spiral_filt_nega{ipatt,t});
                if nansum(temp1(:))~=0
                   temp1_accu_1subject_spiral_neg(:,:,t) = temp1_accu_1subject_spiral_neg(:,:,t) + temp1;
                end
            end
        end
        clear temp1;

        % anticlockwise spirals
        for ipatt = 1:size(spiral_filt_pos,1)
            for t = 1:size(spiral_filt_pos,2)
                temp1 = full(spiral_filt_pos{ipatt,t});
                if nansum(temp1(:))~=0
                   temp1_accu_1subject_spiral_pos(:,:,t) = temp1_accu_1subject_spiral_pos(:,:,t) + temp1;
                end
            end
        end  
        clear temp1;

        t_duration = size(temp1_accu_1subject_spiral_pos,3)-session_duration+1; % last time point to be used for analysis
        % extract task label for each subject
        temp1 = fullTime_allsubject{subject};
        temp1_time = temp1(:,3);
        % find time points for "Right" Tasks (even for hand odd for foot)
        count = find(temp1_time==4); % 4 for Right Hand Tasks
        start_end_time_right_hand = temp1(count,:);
        count = find(temp1_time==3); % 3 for Right Foot Tasks
        start_end_time_right_foot = temp1(count,:);                
        count = find(temp1_time==5); % 5 for tongue cause it is also a right task
        start_end_time_tongue = temp1(count,:);

        % time points for "Left" tasks (even for hand odd for foot)
        count = find(temp1_time==1); % 1 for Left Foot Tasks
        start_end_time_left_foot = temp1(count,:);
        count = find(temp1_time==2); % 2 for Left Hand Tasks
        start_end_time_left_hand = temp1(count,:);                



        % extract task specifc spiral distribution data: "Right" tasks
        % Movement periods for the right hand 
        if size(start_end_time_right_hand,1) == 2
            block_1 = start_end_time_right_hand(1,1):start_end_time_right_hand(1,1)+session_duration;
            block_2 = start_end_time_right_hand(2,1):start_end_time_right_hand(2,1)+session_duration;

            % Task specific spirals grouped into clockwise (-1) and anti-clockwise (+1)
            spiral_distribution_right_hand_Blocks_pos{1,subject} = temp1_accu_1subject_spiral_pos(:,:,block_1);
            spiral_distribution_right_hand_Blocks_pos{2,subject} = temp1_accu_1subject_spiral_pos(:,:,block_2);
            spiral_distribution_right_hand_Blocks_neg{1,subject} = temp1_accu_1subject_spiral_neg(:,:,block_1);
            spiral_distribution_right_hand_Blocks_neg{2,subject} = temp1_accu_1subject_spiral_neg(:,:,block_2);

            %Task specific Spiral centres only
            centres_distribution_right_hand_Blocks_pos{1,subject} = spiral_filt_pos_centre(:,block_1);
            centres_distribution_right_hand_Blocks_pos{2,subject} = spiral_filt_pos_centre(:,block_2);
            centres_distribution_right_hand_Blocks_neg{1,subject} = spiral_filt_nega_centre(:,block_1);
            centres_distribution_right_hand_Blocks_neg{2,subject} = spiral_filt_nega_centre(:,block_2);

            %Task Specific Angular velocity values
            % angVel_profile_right_hand_Blocks_pos{1,subject} = spiral_cav_pos(:,block_1);
            % angVel_profile_right_hand_Blocks_pos{2,subject} = spiral_cav_pos(:,block_2);
            % angVel_profile_right_hand_Blocks_neg{1,subject} = spiral_cav_nega(:,block_1);
            % angVel_profile_right_hand_Blocks_neg{2,subject} = spiral_cav_nega(:,block_2);
            % 
            % Vx_flowmap_right_hand_blocks{1,subject} = Vx_flowmap_norm_phase(:,:,block_1);
            % Vx_flowmap_right_hand_blocks{2,subject} = Vx_flowmap_norm_phase(:,:,block_2);
            % Vy_flowmap_right_hand_blocks{1,subject} = Vy_flowmap_norm_phase(:,:,block_1);
            % Vy_flowmap_right_hand_blocks{2,subject} = Vy_flowmap_norm_phase(:,:,block_2);



            % All detected spiral patterns 
            spiral_patt_right_hand_Blocks_pos{1,subject} = spiral_filt_pos(:,block_1);
            spiral_patt_right_hand_Blocks_pos{2,subject} = spiral_filt_pos(:,block_2);
            spiral_patt_right_hand_Blocks_neg{1,subject} = spiral_filt_nega(:,block_1);
            spiral_patt_right_hand_Blocks_neg{2,subject} = spiral_filt_nega(:,block_2);

        end
       
        
        for block = 1:2
            temp1 = full(spiral_distribution_right_hand_Blocks_neg{block,subject});
            temp2 = full(spiral_distribution_right_hand_Blocks_pos{block,subject});
            if nansum(temp1(:))~=0 && nansum(temp2(:))~=0 %check both conditions to ensure trial equivalence betwen pos and neg
                    no_of_trial = no_of_trial+ 1;
                    neg_right_hand_block_avg(:,:,no_of_trial) = nanmean(temp1(:,:,6:18),3); %spirals of longest duration contribute most to the averages
                    pos_right_hand_block_avg(:,:,no_of_trial) = nanmean(temp2(:,:,6:18),3); 
            end
        end
end  %Subject loop end


    % Get The Sums of the averages For POS and NEG spirals (Heatmap: gives the most common location of long duration spirals)
    neg_right_hand_spiral_heatmap = sum(neg_right_hand_block_avg,3);
    pos_right_hand_spiral_heatmap = sum(pos_right_hand_block_avg,3);


%% Rotational Speed Calculations

for subject = 1
% Initialize arrays to store results for negative and positive vortices
phase_vector_angle_speed_nega = [];
phase_vector_angle_speed_pos = [];

% Limit Analysis Long Lasting spirals >4 frames

% Process negative (clockwise) vortex
for block = [1 2]
    temp1 = spiral_patt_right_hand_Blocks_neg{block,subject};
    for ipatt = 1:size(temp1,1)
        phase_vector_angle_1patt = [];
        % Check if the pattern exists for more than n frames
        if nnz(~cellfun(@isempty,temp1(ipatt,:))) < 4
            continue
        else
        for time = 1:size(temp1,2)
            iSpiral_nega = full(temp1{ipatt,time});
            if isempty(iSpiral_nega) == 1 %Make sure pattern exists at that specific point
                continue
            else %Calcuate the angular speed
            temp2 = centres_distribution_right_hand_Blocks_neg{block,subject};
            iSpiral_nega_centre = temp2{ipatt,time};
            
            % Calculate the filtered phase map
            iSpiral_filt_phase_map = abs(iSpiral_nega) .* smooth_phase_map{subject}(:, :, time);
            
            % Apply a phase filter
            iSpiral_filt_phase_map_pifilt = iSpiral_filt_phase_map;
            iSpiral_filt_phase_map_abs = abs(iSpiral_filt_phase_map);
            iSpiral_filt_phase_map_pifilt(iSpiral_filt_phase_map_abs < 5/8* pi) = 0;

            % Calculate the angle between ~Pi phase points and vortex centers
            idx = find(iSpiral_filt_phase_map_pifilt);
            if isempty(idx)
                phase_vector_angle_1patt = [phase_vector_angle_1patt; nan];
            end
            [row, col] = ind2sub([175, 251], idx);
            phase_vector_u_nega = nanmean(col - iSpiral_nega_centre(1));
            phase_vector_v_nega = nanmean(row - iSpiral_nega_centre(2)); 
            phase_vector_angle_1t1patt = angle(phase_vector_u_nega + 1i * phase_vector_v_nega);
            
            phase_vector_angle_1patt = [phase_vector_angle_1patt; phase_vector_angle_1t1patt];
            end
        end % Time loop End
        % Check if there's enough data for angle speed calculation

        if nansum(phase_vector_angle_1patt(:) ./ phase_vector_angle_1patt(:)) > 1
            % Calculate angle speed and store it
            phase_vector_angle_1patt_diff = anglesubtract(phase_vector_angle_1patt(2:end), phase_vector_angle_1patt(1:end-1));
            phase_vector_angle_speed_nega{block,ipatt} =  phase_vector_angle_1patt_diff;
        end 

        
        for i2 = 1:size(phase_vector_angle_speed_nega{block,ipatt}, 1)
            if phase_vector_angle_speed_nega{block,ipatt}(i2) > pi
                phase_vector_angle_speed_nega{block,ipatt}(i2) = phase_vector_angle_speed_nega{block,ipatt}(i2) - 2 * pi;
            elseif phase_vector_angle_speed_nega{block,ipatt}(i2) < -pi
                phase_vector_angle_speed_nega{block,ipatt}(i2) = phase_vector_angle_speed_nega{block,ipatt}(i2) + 2 * pi;
            end
        end
        end
    end %Pattern Loop End
    try
        non_empty_cells = ~cellfun('isempty',phase_vector_angle_speed_nega(block,:));
    catch ME
        phase_vector_angle_speed_neg{block,1} = [];
        non_empty_cells = 1;
        disp(['subject ' num2str(subject) ' block ' num2str(block) ' is empty'])
    end
        final_phase_vector_angle_speed_nega{block,subject} = phase_vector_angle_speed_nega(block,non_empty_cells);


% Process positive (anti-clockwise) vortex
%empty for subject 4 block 2
    temp1 = spiral_patt_right_hand_Blocks_pos{block,subject};
    for ipatt = 1:size(temp1,1)
        phase_vector_angle_1patt = [];
        % Check if the pattern exists for more than n frames
        if nnz(~cellfun(@isempty,temp1(ipatt,:))) < 4
            continue
        else
        for time = 1:size(temp1,2)
            iSpiral_pos = full(temp1{ipatt,time});
            if isempty(iSpiral_pos) == 1 %Make sure pattern exists at that specific point
                continue
            else %Calcuate the angular speed
            temp2 = centres_distribution_right_hand_Blocks_pos{block,subject};
            iSpiral_pos_centre = temp2{ipatt,time};
            
            % Calculate the filtered phase map
            iSpiral_filt_phase_map = abs(iSpiral_pos) .* smooth_phase_map{subject}(:, :, time);
            
            % Apply a phase filter
            iSpiral_filt_phase_map_pifilt = iSpiral_filt_phase_map;
            iSpiral_filt_phase_map_abs = abs(iSpiral_filt_phase_map);
            iSpiral_filt_phase_map_pifilt(iSpiral_filt_phase_map_abs < 5/8* pi) = 0;

            % Calculate the angle between ~Pi phase points and vortex centers
            idx = find(iSpiral_filt_phase_map_pifilt);
            if isempty(idx)
                phase_vector_angle_1patt = [phase_vector_angle_1patt; nan];
            end
            [row, col] = ind2sub([175, 251], idx);
            phase_vector_u_pos = nanmean(col - iSpiral_pos_centre(1));
            phase_vector_v_pos = nanmean(row - iSpiral_pos_centre(2)); 
            phase_vector_angle_1t1patt = angle(phase_vector_u_pos + 1i * phase_vector_v_pos);
            
            phase_vector_angle_1patt = [phase_vector_angle_1patt; phase_vector_angle_1t1patt];
            end
        end % Time loop End
        % Check if there's enough data for angle speed calculation

        if nansum(phase_vector_angle_1patt(:) ./ phase_vector_angle_1patt(:)) > 1
            % Calculate angle speed and store it
            phase_vector_angle_1patt_diff = anglesubtract(phase_vector_angle_1patt(2:end), phase_vector_angle_1patt(1:end-1));
            phase_vector_angle_speed_pos{block,ipatt} =  phase_vector_angle_1patt_diff;
        end 

        % Regularise to be within [-pi +pi]
        for i2 = 1:size(phase_vector_angle_speed_pos{block,ipatt}, 1)
            if phase_vector_angle_speed_pos{block,ipatt}(i2) > pi
                phase_vector_angle_speed_pos{block,ipatt}(i2) = phase_vector_angle_speed_pos{block,ipatt}(i2) - 2 * pi;
            elseif phase_vector_angle_speed_pos{block,ipatt}(i2) < -pi
                phase_vector_angle_speed_pos{block,ipatt}(i2) = phase_vector_angle_speed_pos{block,ipatt}(i2) + 2 * pi;
            end
        end
        end
    end %Pattern Loop End
    
    try
      non_empty_cells = ~cellfun('isempty',phase_vector_angle_speed_pos(block,:));
    catch ME
      phase_vector_angle_speed_pos{block,1} = [];
      non_empty_cells = 1;
      disp(['subject ' num2str(subject) ' block ' num2str(block) ' is empty'])
    end
      final_phase_vector_angle_speed_pos{block,subject} = phase_vector_angle_speed_pos(block,non_empty_cells);
   
end
end

%    %%  SELECT Most Frequent sprial locations      
%     for irow = 1:size(neg_right_hand_spiral_heatmap,1)
%         for icol = 1:size(neg_right_hand_spiral_heatmap,2)
%             if data > -20   % no of trials contributing to this loc
%                 neg_right_hand_spiral_heatmap(irow,icol) = 0;
%             end
%         end
%     end
% 
%     for irow = 1:size(pos_right_hand_spiral_heatmap,1)
%         for icol = 1:size(pos_right_hand_spiral_heatmap,2)
%             if pos_right_hand_spiral_heatmap(irow,icol) < 20
%                 pos_right_hand_spiral_heatmap(irow,icol) = 0;
%             end
%         end
%     end
%     %%
%     % Calculate Weighted centroids on spiral regions (could be a function)
% 
%      data = nansum(pos_avg(:,:,:),3);
%     for irow = 1:size(data,1)
%         for icol = 1:size(data,2)
%             if data(irow,icol) < 4   % no of trials contributing to this loc
%                 data(irow,icol) = 0;
%             end
%         end
%     end
% 
%     % Create a logical mask for non-zero values
%     nonZeroMask = abs(data) > 0;
% 
%     % Label connected components
%     labelMatrix = bwlabel(nonZeroMask);
% 
%     % Initialize variables to store weighted centroids and total weights
%     weightedCentroids = [];
%     totalWeights = [];
%     disp('run')
%     % Iterate through the regions and extract weighted centroids and total weights
%     for i = 1:max(labelMatrix(:))
%         regionMask = (labelMatrix == i);
%         totalWeight = sum(data(regionMask));
% 
%         % Calculate the weighted centroid for the region
%         [rows, cols] = find(regionMask);
%         weightedCentroidX = sum(rows .* data(regionMask)) / totalWeight;
%         weightedCentroidY = sum(cols .* data(regionMask)) / totalWeight;
% 
%         weightedCentroids = [weightedCentroids; weightedCentroidX, weightedCentroidY];
%         totalWeights = [totalWeights; totalWeight];
%     end
% 
%     % Calculate the final weighted centroid by considering all regions
%     finalWeightedCentroidX = sum(weightedCentroids(:, 1) .* totalWeights) / sum(totalWeights);
%     finalWeightedCentroidY = sum(weightedCentroids(:, 2) .* totalWeights) / sum(totalWeights);
% 
% %% Distance from centroid calc
%     for subject = 1:No_of_Subject
% 
% % NEG SPIRAL DIST
%         for ipatt = 1:size(weightedCentroids,1)
%             avg_spiral_centre = weightedCentroids(ipatt,:);
% 
%             for block = 1:2
%                 temp_dist = [];
%                 temp1 = full(centres_distribution_right_hand_Blocks_neg{block,subject});
%                 emptyRows = all(cellfun(@isempty, temp1),2);
%                 temp1 = temp1(~emptyRows,:); % removes empty rows so no equivalence with ipatt of centre dist
%                 for ipatt2 = 1:size(temp1,1)
%                     centre_timestep = zeros(18,2);
%                     for time = 1:18
%                         if isempty(temp1{ipatt2,time})
% %                           centre_timestep(time,:) = [0 0];
%                             continue
%                         else
%                             centre_timestep(time,:) = temp1{ipatt2,time}; % make sure (x,y)
%                         end
%                     end
%                     if sum(centre_timestep~=0,'all') < 3 %removing 1 frame spirals  
%                         continue
%                     else
%                         distance_timestep = zeros(18,1);
%                         for time = 1:18
%                             if centre_timestep(time,1) == 0
%                                 continue
%                             else
%                             distance_timestep(time) = sqrt((centre_timestep(time,2) - avg_spiral_centre(2))^2+(centre_timestep(time,1) - avg_spiral_centre(1))^2);
%                             end
%                         end
% 
%                     %eff_spiral_centre = nanmean(centre_timestep,1);
%                     % try average or closest distance
% %                     index = find(centre_timestep,1,'last')-18;   %18 Window duration
% %                     eff_spiral_centre = centre_timestep(index,:);
% 
% 
%                     % Calculate the distances from each point in A to leftbottom
%                  %  temp_dist = [temp_dist sqrt((eff_spiral_centre(2) - avg_spiral_centre(2))^2+(eff_spiral_centre(1) - avg_spiral_centre(1))^2)];                        
%                     temp_dist = [temp_dist  min(distance_timestep(distance_timestep>0))];                        
%                     % Find the index of the minimum distance
%                     [~, min_index] = min(temp_dist);
%                     end
%                 end
%                 distances_neg(subject,block,ipatt)=(min(temp_dist));
%             end
%         end
        
%         % POS SPIRAL DIST
%         for ipatt = 1:size(weightedCentroids,1)
%             avg_spiral_centre = weightedCentroids(ipatt,:); %(y,x)
%             
%             for block = 1:2
%                 temp_dist = [];
%                 temp1 = full(centres_distribution_right_hand_Blocks_pos{block,subject});
%                 emptyRows = all(cellfun(@isempty, temp1),2);
%                 temp1 = temp1(~emptyRows,:); % removes empty rows so no equivalence with ipatt of centre dist
%                 for ipatt2 = 1:size(temp1,1)
%                     centre_timestep = zeros(18,2);
%                     for time = 1:18
%                         if isempty(temp1{ipatt2,time})
% %                           centre_timestep(time,:) = [0 0];
%                             continue
%                         else
%                             centre_timestep(time,:) = temp1{ipatt2,time}; % make sure (x,y)
%                             
% %                             hold on
% %                             x =  full(spiral_distribution_right_hand_Blocks_pos{block,subject});
% %                             imagesc(x(:,:,time))
% %                             scatter(centre_timestep(time,:),'white')
% %                             set(gca,'ydir','normal')
% %                             colormap jet
% %                             pause(20)
%                             
%                         end
%                     end
%                     if sum(centre_timestep~=0,'all') < 3 %removing 1 frame spirals  
%                         continue
%                     else
%                         
%                         distance_timestep = zeros(18,1);
%                         for time = 1:18
%                             if centre_timestep(time,1) == 0
%                                 continue
%                             else
%                             distance_timestep(time) = sqrt((centre_timestep(time,1) - avg_spiral_centre(2))^2+(centre_timestep(time,2) - avg_spiral_centre(1))^2);
%                             end
%                         end
% 
%                     %eff_spiral_centre = nanmean(centre_timestep,1);
%                     % try average or closest distance
% %                     index = find(centre_timestep,1,'last')-18;   %18 Window duration
% %                     eff_spiral_centre = centre_timestep(index,:);
% 
%                 
%                     % Calculate the distances from each point in A to leftbottom
%                  %  temp_dist = [temp_dist sqrt((eff_spiral_centre(1) - avg_spiral_centre(2))^2+(eff_spiral_centre(2) - avg_spiral_centre(1))^2)];                        
%                     temp_dist = [temp_dist  min(distance_timestep(distance_timestep>0))];                        
%                     % Find the index of the minimum distance
%                     [~, min_index] = min(temp_dist);
%                     end
%                 end
%                 distances_pos(subject,block,ipatt)=(min(temp_dist));
%             end
%         end
% end


% %% Visualisation of Heat Map
%         figure()
%         hold on
% 
%         %imagesc(pos_right_hand_block_avg(:,:,1)) 
%         %imagesc(pos_right_hand_spiral_heatmap+neg_right_hand_spiral_heatmap) 
%         imagesc(data)
%         %plot(weightedCentroids(:,2),weightedCentroids(:,1),'white')
%         set(gca,'ydir','normal')
%         colormap jet
% 
%         for parcellation_ID = 1:22 % Numbers correspond to certain areas of brain
%             parcellation_template_1par = parcellation_template;
%             parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
% 
%             parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
% 
%             B = bwboundaries(parcellation_template_1par(:,:,1),'noholes');
%             for k = 1:length(B)
%                 boundary = B{k};
%                 plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])
% 
%             end
%         end 
%         hold off
%         colorbar
% %%
% no_of_trial = 0;
% for subject = [5 8 10 11 13 18 21 24 35 41 42 44 45 64 67 77 82 84 95 99 102 103 114 136]
% 
%         for block = 1:2
%             temp1 = full(spiral_distribution_right_hand_Blocks_neg{block,subject});
%             temp2 = full(spiral_distribution_right_hand_Blocks_pos{block,subject});
%             if nansum(temp1(:))~=0 && nansum(temp2(:))~=0 %check both conditions to ensure trial equivalence betwen pos and neg
%                     no_of_trial = no_of_trial+ 1;
%                     neg_avg(:,:,no_of_trial) = nanmean(temp1(:,:,6:18),3); %spirals of longest duration contribute most to the averages
%                     pos_avg(:,:,no_of_trial) = nanmean(temp2(:,:,6:18),3); 
%             end
%         end
% end
% 
% first = nansum(neg_avg(:,:,25:48),3);
% %nansum(pos_avg,3);
% 
% figure()
% hold on
% imagesc(first)
% set(gca,'ydir','normal')
% colormap jet   
% colorbar()
% for parcellation_ID = 1:22 % Numbers correspond to certain areas of brain
%     parcellation_template_1par = parcellation_template;
%     parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
% 
%     parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
% 
%     B = bwboundaries(parcellation_template_1par(:,:,1),'noholes');
%     for k = 1:length(B)
%         boundary = B{k};
%         plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])
% 
%     end
% end 
% 
% 
 %%

           folder_name = [main_folder,'/Sample Data/Motor Task/Analysis'];
           cd(folder_name)
           filename = ['rotational_speeds','.mat'];
           save([folder_name,'/',filename]);     

end

