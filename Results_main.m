%% Analysis Procedure for main results
clear all
% profile on

restoredefaultpath
% cd '/headnode2/akar5239/BrainVortexToolbox-main/'
cd 'E:\resarch_data\fMRI\BrainVortexToolbox-main'

addpath(genpath([pwd]))
main_folder = pwd;

flagSur = 1;%  0 for real data, 1 to generate surrogate data

% a range of parameters avaiable for different dataset, but for demonstration
% purpose, only use the parameters provided
No_of_Subject = 22; % number of subjects used for analysis, randomly selected from HCP database (S1200)
flagRest = 0; %  resting data => 1 , task data => 0 

hemisphere = 1; % 1 for left hemisphere, 2 for right hemisphere, keep this as 1 for demonstration purpose

% flagTask: 
% 1 = language task, original 100 subjects; 
% 2 = language task, additional 100 subjects;
% 3 = working memory task;
% 4 = Motor Task
flagTask = 4;  

% 0 = unsmoothed, raw data; 1 = temporally smoothed (bandpass filtered)
% data; 2 = spatiotemporally smoothed (bandpass filtered) data
flagSmooth = 1;


% Preprocessing of raw fMRI data from HCP 
for flagSur = 0:1; %  0 for real data, 1 to generate surrogate data, both are needed for further analysis
    for hemisphere = 1:1; % 1 for left hemisphere, 2 for right hemisphere
        % flagTask: 1 = language task, original 100 subjects; 2 = language task, additional 100 subjects
        % 3 = working memory task
        for flagTask = 4:4
        % 0 = unsmoothed, raw data; 1 = temporally smoothed (bandpass filtered)
        % data; 2 = spatiotemporally smoothed (bandpass filtered) data
            for flagSmooth = 1:2; 
                for subject = 1:No_of_Subject
                    fprintf("Preprocessing for sub" + subject + "\n")
                    [sigBPass] = Preprocessing_main(subject,hemisphere,flagSur,flagRest,flagTask,flagSmooth);
                end     
            end
        end
    end
end


%% Full-size Vortex detection followed by statistical testing against sprials detected in the null model (surrogate data)
% only sprials that contain similarity index higher than 95th percentile of
% the null model (surrogate data) are seleceted for further analysis
% currentsub = 26;
flagRest = 0;
for flagTask = 4:4
    fprintf('flagTask=%d\n',flagTask)
    for hemisphere = 1:1
        for subject = 1:No_of_Subject
            fprintf("Detecting spirals for sub" + subject + "\n")
            if flagRest == 0
                if flagTask == 1
                    if hemisphere == 1
                       % load spatiotemporally bandpass filtered (smoothed) real data
                       folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        DataIn_smooth = permute(DataOut(:,1:175,:,1:315),[2,3,4,1]);;;
                       DataIn_smooth = DataOut(1:175,:,1:315);
                       % load temporally bandpass filtered (unsmoothed) real data
                       folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_unsmooth = DataOut(1:175,:,1:315);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_sur_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;   
                    elseif hemisphere == 2
                       % load spatiotemporally bandpass filtered (smoothed) real data
                       folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_language_task_orig100_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        DataIn_smooth = permute(DataOut(:,1:175,:,1:315),[2,3,4,1]);
                       DataIn_smooth = DataOut(1:175,:,1:315);
                       % load temporally bandpass filtered (unsmoothed) real data
                       folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_language_task_orig100_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_unsmooth = DataOut(1:175,:,1:315);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_language_task_orig100_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Language Task Original 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_sur_language_task_orig100_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;                   
                    end
                elseif flagTask == 2
                    if hemisphere == 1
                       % load spatiotemporally bandpass filtered (smoothed) real data
                       folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        DataIn_smooth = permute(DataOut(:,1:175,:,1:315),[2,3,4,1]);;
                       DataIn_smooth = DataOut(1:175,:,1:315);
                       % load temporally bandpass filtered (unsmoothed) real data
                       folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_unsmooth = DataOut(1:175,:,1:315);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_sur_language_task_orig100_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;   
                    elseif hemisphere == 2
                       % load spatiotemporally bandpass filtered (smoothed) real data
                       folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_language_task_orig100_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        DataIn_smooth = permute(DataOut(:,1:175,:,1:315),[2,3,4,1]);;
                       DataIn_smooth = DataOut(1:175,:,1:315);
                       % load temporally bandpass filtered (unsmoothed) real data
                       folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_language_task_orig100_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_unsmooth = DataOut(1:175,:,1:315);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_language_task_orig100_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Language Task Additional 100 sub/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_sur_language_task_orig100_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;                   
                    end 
                elseif flagTask == 3
                    if hemisphere == 1
                       % load spatiotemporally bandpass filtered (smoothed) real data
                       folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_WM_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        DataIn_smooth = permute(DataOut(:,1:175,:,:),[2,3,4,1]);
                       DataIn_smooth = DataOut(1:175,:,:);
                       % load temporally bandpass filtered (unsmoothed) real data
                       folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_WM_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_unsmooth = DataOut(1:175,:,:);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_WM_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_sur_WM_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;   
                    elseif hemisphere == 2
                       % load spatiotemporally bandpass filtered (smoothed) real data
                       folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_WM_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        DataIn_smooth = permute(DataOut(:,1:175,:,:),[2,3,4,1]);;
                       DataIn_smooth = DataOut(1:175,:,1:315);
                       % load temporally bandpass filtered (unsmoothed) real data
                       folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_WM_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_unsmooth = DataOut(1:175,:,1:315);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_WM_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Working Memory Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_sur_WM_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;   
                    end
                elseif flagTask == 4
                    if hemisphere == 1
                       % load spatiotemporally bandpass filtered (smoothed) real data
                       folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        
                    %    DataIn_smooth = DataOut(1:175,:,1:283);
                       DataIn_smooth = DataOut(1:175,:,1:217); % 217 for LSD data
                       % load temporally bandpass filtered (unsmoothed) real data
                       folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_motor_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                    %    DataIn_unsmooth = DataOut(1:175,:,1:283);
                       DataIn_unsmooth = DataOut(1:175,:,1:217);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_motor_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_sur_motor_task_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;   
                    elseif hemisphere == 2
                       % load spatiotemporally bandpass filtered (smoothed) real data
                       folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_motor_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        DataIn_smooth = permute(DataOut(:,1:175,:,1:315),[2,3,4,1]);
                       DataIn_smooth = DataOut(1:175,:,1:315);
                       % load temporally bandpass filtered (unsmoothed) real data
                       folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_motor_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_unsmooth = DataOut(1:175,:,1:315);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_motor_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Motor Task/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_sur_motor_task_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;                   
                    end
                end
            elseif flagRest == 1
                    if hemisphere == 1
                       % load spatiotemporally bandpass filtered (smoothed) real data
                       folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_resting_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        DataIn_smooth = permute(DataOut(:,1:175,:,1:1199),[2,3,4,1]);;;
                       DataIn_smooth = DataOut(1:175,:,1:1199);
                       % load temporally bandpass filtered (unsmoothed) real data
                       folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_resting_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_unsmooth = DataOut(1:175,:,1:1199);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_resting_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_sur_resting_LEFT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;   
                    elseif hemisphere == 2
                       % load spatiotemporally bandpass filtered (smoothed) real data
                       folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_resting_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
%                        DataIn_smooth = permute(DataOut(:,1:175,:,1:1199),[2,3,4,1]);;;
                       DataIn_smooth = DataOut(1:175,:,1:1199);
                       % load temporally bandpass filtered (unsmoothed) real data
                       folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_resting_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_unsmooth = DataOut(1:175,:,1:1199);
                       % load spatiotemporally bandpass filtered (smoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_spatiotemporalbandpass_data_sur_resting_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_smooth = DataOut_smooth;         
                       % load temporally bandpass filtered (unsmoothed) surrogate data
                       folder_name = [main_folder,'/Sample Data/Resting/Preprocessed Data'];
                       cd(folder_name)
                       filename = ['Preprocessed_temporalbandpass_data_sur_resting_RIGHT_sub',num2str(subject),'.mat'];
                       load(filename)
                       DataIn_sur_unsmooth = DataOut_unsmooth;                     
                    end
            end
            [temp1_compatibility_ratio_pos_nega_real_avg,temp1_compatibility_ratio_pos_nega_sur_avg] = spiral_detection_surfilt(subject,DataIn_smooth,DataIn_unsmooth, DataIn_sur_smooth, DataIn_sur_unsmooth,main_folder,flagRest,flagSmooth,flagTask,hemisphere)
        end
    end
end

%% Spiral distribution zscore map: Fig .2d & radii & lifetime & propagation speed 

% [spiral_template_timeavg_accu_stdnorm, spiral_radius_accu_avg, spiral_duration_accu_avg, spiral_count_timeavg, spiral_transverse_speed_accu_avg] = spiral_distribution_zscore_map_speed_duration_radius_count(flagRest,flagSur,flagSmooth,flagTask,hemisphere,main_folder,No_of_Subject);

%% fMRI amplitude vs. distance from singularity: Fig 2e

% [distance_from_centre_amplitude_data_pos_accu] = distance_vs_amplitude(No_of_Subject,flagRest,flagTask,flagSur,flagSmooth,hemisphere,main_folder);

%% sprial interaction statistics: Fig 4a

% [count_repulsion_percent_avg,count_partialannihilation_percent_avg,count_fullannihilation_percent_avg] = spiral_interaction_statistics(flagRest,flagSur,hemisphere,main_folder,No_of_Subject,flagTask);

%% extract task-label (language task and working memory task)

%%%%%%%%%%%%%%%%%%%%%%%[task_type_name] = TaskLabel_Extract(flagRest,main_folder,No_of_Subject,flagTask);

%% Task-specific trial-averaged spiral distribution map: Fig 5 (language task)
%%%%%%%%%%%%%%%%%%%%%%%flagTask = 4; % keep this as 1 for this demonstration

%%%%%%%%%%%%%%%%%%%%%%%[spiral_distribution_math_listen_mid_avg,spiral_distribution_math_answer_end_avg,spiral_distribution_story_listen_mid_avg,spiral_distribution_story_answer_end_avg] = task_specific_spiral_distribution(flagRest,flagSur,hemisphere,main_folder,No_of_Subject,flagTask);

%% Task-specific trial-averaged spiral distribution, contrast map: Fig 5 (language task)
% and contrast significance distributions in 7 functional networks: Fig 6a
%flagTask = 1; % keep this as 1 for this demonstration

%[p_value_negative_log_math_story_listen,p_value_negative_log_math_story_answer] = spiral_contrast_significance(flagRest,flagSur,hemisphere,main_folder,No_of_Subject,flagTask)

%% sprial centre-based classifier: Fig 6b and c (language task)
%flagTask = 1; % keep this as 1 for this demonstration

%[classification_accuracy_avg, classification_accuracy_stderr, classification_accuracy] = spiral_classifer_language_task(flagRest,flagSur,hemisphere,main_folder,No_of_Subject,flagTask);

%% task-evoked trial-averaged unfiltered fMRI signal: Fig 7a-c (language task)
%listen_or_answer = 1; % 1 for listening tasks; 2 for answering tasks
%flagTask = 1; % keep this as 1 for this demonstration

%[temp1_phase,Vx_flowmap_norm_phase,Vy_flowmap_norm_phase] = task_evoked_unfiltered_fMRI_signal(flagSur,hemisphere,main_folder,No_of_Subject,flagTask,listen_or_answer);

%% PCA analysis: Fig 7d (language task)
%listen_or_answer = 1; % 1 for listening tasks; 2 for answering tasks
%motor_or_PCC = 1; % 1 for motor (M1-PMd); 2 for PCC
%flagTask = 1; % keep this as 1 for this demonstration

%[score,latent] = PCA(flagSur,hemisphere,main_folder,No_of_Subject,flagTask,listen_or_answer,motor_or_PCC);

%% Region of coordination (ROC): Fig 8b & 8d (language task)
%flagTask = 1; % keep this as 1 for this demonstration

%[region_of_coordination] = region_of_coordination_ROC(flagSur,hemisphere,main_folder,No_of_Subject,flagTask);

%% local phase vector field based classifier: Fig 8c (language task)
%flagTask = 1; % keep this as 1 for this demonstration

%[classification_accuracy_avg,classification_accuracy_stderr] = PhaseVectorAngle_local_classifier(flagSur,hemisphere,main_folder,No_of_Subject,flagTask);


% profile off
% profile viewer
% 
% profile -history
% P = profile('info');
% save('profile_detecting_spirals.mat', 'P');