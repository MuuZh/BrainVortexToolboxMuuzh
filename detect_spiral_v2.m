function  [temp1_compatibility_ratio_pos_nega_real_avg,temp1_compatibility_ratio_pos_nega_sur_avg]=detect_spiral_v2(subject_index)

%% Analysis Procedure for main results
    warning('off','all')
    cd '/headnode2/mtan8991/BrainVortexToolboxMuuzh-main/'
   %  cd 'E:\resarch_data\fMRI\BrainVortexToolboxMuuzh-main'
    addpath(genpath([pwd]))
    main_folder = pwd;
    
    flagSur = 1;%  0 for real data, 1 to generate surrogate data
    
    % a range of parameters avaiable for different dataset, but for demonstration
    % purpose, only use the parameters provided
    flagRest = 0; %  resting data => 1 , task data => 0 
    
    hemisphere = 2; % 1 for left hemisphere, 2 for right hemisphere, keep this as 1 for demonstration purpose
    
    % flagTask: 
    % 1 = language task, original 100 subjects; 
    % 2 = language task, additional 100 subjects;
    % 3 = working memory task;
    % 4 = Motor Task
    flagTask = 4;  
    
    % 0 = unsmoothed, raw data; 1 = temporally smoothed (bandpass filtered)
    % data; 2 = spatiotemporally smoothed (bandpass filtered) data
    flagSmooth = 1;

   task_folder_path = '/data/DMT_data/PCB_group';
   working_directory_path = [main_folder, task_folder_path];
   raw_data_path = [working_directory_path,'/Raw Data'];
   data_pos_path = [working_directory_path,'/Data Pos'];
   preprocessed_path = [working_directory_path,'/Preprocessed Data'];
   spiral_detected_path = [working_directory_path,'/Spiral Detected'];


   raw_data_dir = dir(raw_data_path);
   data_pos_dir = dir(data_pos_path);
   preprocessed_dir = dir(preprocessed_path);
%% Full-size Vortex detection followed by statistical testing against sprials detected in the null model (surrogate data)
% only sprials that contain similarity index higher than 95th percentile of
% the null model (surrogate data) are seleceted for further analysis
flagRest = 0;
for flagTask = 4
    for hemisphere = 1
        for subject = subject_index 
            if flagRest == 0
                if flagTask == 4
                  cd(preprocessed_path)
                    if hemisphere == 1
                       % load spatiotemporally bandpass filtered (smoothed) real data
                       filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        
                       DataIn_smooth = DataOut(1:175,:,1:217);
                       % load temporally bandpass filtered (unsmoothed) real data
                       folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];

                       filename = ['Preprocessed_temporalbandpass_data_motor_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_unsmooth = DataOut(1:175,:,1:217);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data

                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_motor_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data

                       filename = ['Preprocessed_temporalbandpass_data_sur_motor_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;   
                    elseif hemisphere == 2
                       % load spatiotemporally bandpass filtered (smoothed) real data

                       filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        DataIn_smooth = permute(DataOut(:,1:175,:,1:315),[2,3,4,1]);
                       DataIn_smooth = DataOut(1:175,:,1:217);
                       % load temporally bandpass filtered (unsmoothed) real data

                       filename = ['Preprocessed_temporalbandpass_data_motor_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_unsmooth = DataOut(1:175,:,1:217);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data

                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_motor_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data

                       filename = ['Preprocessed_temporalbandpass_data_sur_motor_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;                   
                    end
                end
            elseif flagRest == 1
                    1;
            end
            [temp1_compatibility_ratio_pos_nega_real_avg,temp1_compatibility_ratio_pos_nega_sur_avg] = spiral_detection_surfilt(subject,DataIn_smooth,DataIn_unsmooth, DataIn_sur_smooth, DataIn_sur_unsmooth,main_folder,flagRest,flagSmooth,flagTask,hemisphere);
        end
    end
end
end
