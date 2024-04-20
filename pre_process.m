function [sigBPass] = pre_process(subject_index)
    %% Analysis Procedure for main results
    warning('off','all')
    cd '/headnode2/mtan8991/BrainVortexToolbox-main/'
    % cd 'E:\resarch_data\fMRI\BrainVortexToolbox-main'
    % cd 'E:/resarch_data/fMRI/'
    addpath(genpath([pwd]))
    main_folder = pwd;
    
    flagSur = 1;%  0 for real data, 1 to generate surrogate data
    
    % a range of parameters avaiable for different dataset, but for demonstration
    % purpose, only use the parameters provided
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
    
    
    %% Preprocessing of raw fMRI data from HCP 
    for flagSur = 0:1 %  0 for real data, 1 to generate surrogate data, both are needed for further analysis
        for hemisphere = 1 % 1 for left hemisphere, 2 for right hemisphere
            % flagTask: 1 = language task, original 100 subjects; 2 = language task, additional 100 subjects
            % 3 = working memory task
            for flagTask = 4
            % 0 = unsmoothed, raw data; 1 = temporally smoothed (bandpass filtered)
            % data; 2 = spatiotemporally smoothed (bandpass filtered) data
                for flagSmooth = 1:2
                    for subject = subject_index 
                        fprintf("Preprocessing for sub" + subject + "\n")
                        [sigBPass] = Preprocessing_main(subject,hemisphere,flagSur,flagRest,flagTask,flagSmooth);
                    end     
                end
            end
        end
    end


end
