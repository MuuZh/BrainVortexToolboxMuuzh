function [subject] = new_task_specific_spiral(subject)
restoredefaultpath
cd '/headnode2/akar5239/BrainVortexToolbox-main/'
addpath(genpath([pwd]))
main_folder = pwd;

flagSur = 0;%  0 for real data, 1 to generate surrogate data

% a range of parameters avaiable for different dataset, but for demonstration
% purpose, only use the parameters provided
No_of_Subject = 1; % number of subjects used for analysis, randomly selected from HCP database (S1200)
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


%% load task label (language task)
disp(['loading task label...'])
if flagTask == 4 % Motor Task
    % load task label of each subject`
    foldername = [main_folder,'/Sample Data/Motor Task/TaskLabel'];
    cd(foldername)
    name = dir(pwd) ;
    file_name2 = ['MotorTaskLabelAllSubject.mat'];
    load (file_name2);
    for isubject = 1:150
        fullTime_allsubject{isubject} = TaskLabel_AllSubject_motor{isubject};
    end
end



%% Motor Task
if flagTask == 4 % Motor Task
    % define parameters
    session_duration = 17;   % (currently this is the number of frames 12 seconds/0.72 )
    disp(['extracting task specific spiral centre distributions...'])


    foldername = [main_folder,'/Sample Data/Motor Task/Spiral Detected'];
    cd(foldername)
    if hemisphere == 1
        filename = ['Spiral_detected_surfilt_motor_task_LEFT_sub',num2str(subject),'.mat'];
        load(filename)
    elseif hemisphere == 2
        filename = ['Spiral_detected_surfilt_motor_task_RIGHT_sub',num2str(subject),'.mat'];
        load(filename)
    end

    if flagSur == 0
        spiral_filt_nega = spiral_filt_nega_real;
        spiral_filt_pos = spiral_filt_pos_real;

        % spiral_filt_nega = spiral_filt_nega_real_95perc_extend;
        % spiral_filt_pos = spiral_filt_pos_real_95perc_extend;
    elseif flagSur == 1
        spiral_filt_nega = spiral_filt_nega_sur;
        spiral_filt_pos = spiral_filt_pos_sur;
    end

    if flagRest == 0
        temp1_accu_1subject_spiral = zeros(175,251,283);
    elseif flagRest == 1
        temp1_accu_1subject_spiral = zeros(175,251,1199);
    end


    % extract task label for each subject
    temp1 = fullTime_allsubject{subject};
    temp1_time = temp1(:,3);
    % find time points for "Right" Tasks (even for hand odd for foot)
    count = find(temp1_time==4); % 4 for Right Hand Tasks
    start_end_time_right_hand = temp1(count,:);
    rh_cues =  [start_end_time_right_hand(1,1)-4; start_end_time_right_hand(2,1)-4];        %THIS NEEDS TO BE FIXED
    count = find(temp1_time==3); % 3 for Right Foot Tasks
    start_end_time_right_foot = temp1(count,:);
    rf_cues =  [start_end_time_right_foot(1,1)-4; start_end_time_right_foot(2,1)-4];
    count = find(temp1_time==5); % 5 for tongue cause it is also a right task
    start_end_time_tongue = temp1(count,:);
    t_cues =  [start_end_time_tongue(1,1)-4; start_end_time_tongue(2,1)-4];

    % time points for "Left" tasks (even for hand odd for foot)
    count = find(temp1_time==1); % 1 for Left Foot Tasks
    start_end_time_left_foot = temp1(count,:);
    lf_cues =  [start_end_time_left_foot(1,1)-4; start_end_time_left_foot(2,1)-4];
    count = find(temp1_time==2); % 2 for Left Hand Tasks
    start_end_time_left_hand = temp1(count,:);
    lh_cues =  [start_end_time_left_hand(1,1)-4; start_end_time_left_hand(2,1)-4];





    % extract task specifc spiral distribution data: "Right" tasks
    % Movement periods for the right hand
    if size(start_end_time_right_hand,1) == 2
        block_1 = start_end_time_right_hand(1,1):start_end_time_right_hand(1,1)+session_duration;
        block_2 = start_end_time_right_hand(2,1):start_end_time_right_hand(2,1)+session_duration;


        [spontaneuos_vortices, spontaneuous_vortex_longest_duration, modulated_vortices, modulated_vortices_longest_duration] = task_onset_id(spiral_filt_pos,spiral_filt_pos_centreONLY_real,rh_cues, [block_1; block_2]);

        rh_block1_pos_longest_spontaneous = spontaneuous_vortex_longest_duration(1);
        rh_block2_pos_longest_spontaneous = spontaneuous_vortex_longest_duration(2);

        rh_block1_pos_longest_modulated = modulated_vortices_longest_duration(1);
        rh_block2_pos_longest_modulated = modulated_vortices_longest_duration(2);


        rh_block1_pos_spontaneous_vortex = spontaneuos_vortices{1};
        rh_block2_pos_spontaneous_vortex = spontaneuos_vortices{2};

        rh_block1_pos_modulated_vortex = modulated_vortices{1};
        rh_block2_pos_modulated_vortex = modulated_vortices{2};

        [spontaneuos_vortices, spontaneuous_vortex_longest_duration, modulated_vortices, modulated_vortices_longest_duration] = task_onset_id(spiral_filt_nega,spiral_filt_nega_centreONLY_real, rh_cues, [block_1; block_2]);


        rh_block1_neg_longest_spontaneous = spontaneuous_vortex_longest_duration(1);
        rh_block2_neg_longest_spontaneous = spontaneuous_vortex_longest_duration(2);

        rh_block1_neg_longest_modulated = modulated_vortices_longest_duration(1);
        rh_block2_neg_longest_modulated = modulated_vortices_longest_duration(2);


        rh_block1_nega_spontaneous_vortex = spontaneuos_vortices{1};
        rh_block2_nega_spontaneous_vortex = spontaneuos_vortices{2};

        rh_block1_nega_modulated_vortex = modulated_vortices{1};
        rh_block2_nega_modulated_vortex = modulated_vortices{2};
    end

    % Movement periods for the right foot
    if size(start_end_time_right_foot,1) == 2
        block_1 = start_end_time_right_foot(1,1):start_end_time_right_foot(1,1)+session_duration;
        block_2 = start_end_time_right_foot(2,1):start_end_time_right_foot(2,1)+session_duration;

        [spontaneuos_vortices, spontaneuous_vortex_longest_duration, modulated_vortices, modulated_vortices_longest_duration] = task_onset_id(spiral_filt_pos, spiral_filt_pos_centreONLY_real ,rf_cues, [block_1; block_2]);

        rf_block1_pos_longest_spontaneous = spontaneuous_vortex_longest_duration(1);
        rf_block2_pos_longest_spontaneous = spontaneuous_vortex_longest_duration(2);

        rf_block1_pos_longest_modulated = modulated_vortices_longest_duration(1);
        rf_block2_pos_longest_modulated = modulated_vortices_longest_duration(2);


        rf_block1_pos_spontaneous_vortex = spontaneuos_vortices{1};
        rf_block2_pos_spontaneous_vortex = spontaneuos_vortices{2};

        rf_block1_pos_modulated_vortex = modulated_vortices{1};
        rf_block2_pos_modulated_vortex = modulated_vortices{2};

        [spontaneuos_vortices, spontaneuous_vortex_longest_duration, modulated_vortices, modulated_vortices_longest_duration] = task_onset_id(spiral_filt_nega,spiral_filt_nega_centreONLY_real, rf_cues, [block_1; block_2]);


        rf_block1_neg_longest_spontaneous = spontaneuous_vortex_longest_duration(1);
        rf_block2_neg_longest_spontaneous = spontaneuous_vortex_longest_duration(2);

        rf_block1_neg_longest_modulated = modulated_vortices_longest_duration(1);
        rf_block2_neg_longest_modulated = modulated_vortices_longest_duration(2);


        rf_block1_nega_spontaneous_vortex = spontaneuos_vortices{1};
        rf_block2_nega_spontaneous_vortex = spontaneuos_vortices{2};

        rf_block1_nega_modulated_vortex = modulated_vortices{1};
        rf_block2_nega_modulated_vortex = modulated_vortices{2};

    end

    % Movement periods for the tongue
    if size(start_end_time_tongue,1) == 2
        block_1 = start_end_time_tongue(1,1):start_end_time_tongue(1,1)+session_duration;
        block_2 = start_end_time_tongue(2,1):start_end_time_tongue(2,1)+session_duration;

        [spontaneuos_vortices, spontaneuous_vortex_longest_duration, modulated_vortices, modulated_vortices_longest_duration] = task_onset_id(spiral_filt_pos,spiral_filt_pos_centreONLY_real ,t_cues, [block_1; block_2]);

        t_block1_pos_longest_spontaneous = spontaneuous_vortex_longest_duration(1);
        t_block2_pos_longest_spontaneous = spontaneuous_vortex_longest_duration(2);

        t_block1_pos_longest_modulated = modulated_vortices_longest_duration(1);
        t_block2_pos_longest_modulated = modulated_vortices_longest_duration(2);


        t_block1_pos_spontaneous_vortex = spontaneuos_vortices{1};
        t_block2_pos_spontaneous_vortex = spontaneuos_vortices{2};

        t_block1_pos_modulated_vortex = modulated_vortices{1};
        t_block2_pos_modulated_vortex = modulated_vortices{2};

        [spontaneuos_vortices, spontaneuous_vortex_longest_duration, modulated_vortices, modulated_vortices_longest_duration] = task_onset_id(spiral_filt_nega, spiral_filt_nega_centreONLY_real,t_cues, [block_1; block_2]);


        t_block1_neg_longest_spontaneous = spontaneuous_vortex_longest_duration(1);
        t_block2_neg_longest_spontaneous = spontaneuous_vortex_longest_duration(2);

        t_block1_neg_longest_modulated = modulated_vortices_longest_duration(1);
        t_block2_neg_longest_modulated = modulated_vortices_longest_duration(2);


        t_block1_nega_spontaneous_vortex = spontaneuos_vortices{1};
        t_block2_nega_spontaneous_vortex = spontaneuos_vortices{2};

        t_block1_nega_modulated_vortex = modulated_vortices{1};
        t_block2_nega_modulated_vortex = modulated_vortices{2};
    end

    % extract task specifc spiral distribution data: "Left" tasks
    % Movement Periods for the Left Foot
    if size(start_end_time_left_foot,1) == 2
        block_1 = start_end_time_left_foot(1,1):start_end_time_left_foot(1,1)+session_duration;
        block_2 = start_end_time_left_foot(2,1):start_end_time_left_foot(2,1)+session_duration;

        [spontaneuos_vortices, spontaneuous_vortex_longest_duration, modulated_vortices, modulated_vortices_longest_duration] = task_onset_id(spiral_filt_pos, spiral_filt_pos_centreONLY_real,lf_cues, [block_1; block_2]);

        lf_block1_pos_longest_spontaneous = spontaneuous_vortex_longest_duration(1);
        lf_block2_pos_longest_spontaneous = spontaneuous_vortex_longest_duration(2);

        lf_block1_pos_longest_modulated = modulated_vortices_longest_duration(1);
        lf_block2_pos_longest_modulated = modulated_vortices_longest_duration(2);


        lf_block1_pos_spontaneous_vortex = spontaneuos_vortices{1};
        lf_block2_pos_spontaneous_vortex = spontaneuos_vortices{2};

        lf_block1_pos_modulated_vortex = modulated_vortices{1};
        lf_block2_pos_modulated_vortex = modulated_vortices{2};

        [spontaneuos_vortices, spontaneuous_vortex_longest_duration, modulated_vortices, modulated_vortices_longest_duration] = task_onset_id(spiral_filt_nega,spiral_filt_nega_centreONLY_real ,lf_cues, [block_1; block_2]);

        lf_block1_neg_longest_spontaneous = spontaneuous_vortex_longest_duration(1);
        lf_block2_neg_longest_spontaneous = spontaneuous_vortex_longest_duration(2);

        lf_block1_neg_longest_modulated = modulated_vortices_longest_duration(1);
        lf_block2_neg_longest_modulated = modulated_vortices_longest_duration(2);


        lf_block1_nega_spontaneous_vortex = spontaneuos_vortices{1};
        lf_block2_nega_spontaneous_vortex = spontaneuos_vortices{2};

        lf_block1_nega_modulated_vortex = modulated_vortices{1};
        lf_block2_nega_modulated_vortex = modulated_vortices{2};

    end


    if size(start_end_time_left_hand,1) == 2
        block_1 = start_end_time_left_hand(1,1):start_end_time_left_hand(1,1)+session_duration;
        block_2 = start_end_time_left_hand(2,1):start_end_time_left_hand(2,1)+session_duration;

        [spontaneuos_vortices, spontaneuous_vortex_longest_duration, modulated_vortices, modulated_vortices_longest_duration] = task_onset_id(spiral_filt_pos, spiral_filt_pos_centreONLY_real,lh_cues, [block_1; block_2]);

        lh_block1_pos_longest_spontaneous = spontaneuous_vortex_longest_duration(1);
        lh_block2_pos_longest_spontaneous = spontaneuous_vortex_longest_duration(2);

        lh_block1_pos_longest_modulated = modulated_vortices_longest_duration(1);
        lh_block2_pos_longest_modulated = modulated_vortices_longest_duration(2);


        lh_block1_pos_spontaneous_vortex = spontaneuos_vortices{1};
        lh_block2_pos_spontaneous_vortex = spontaneuos_vortices{2};

        lh_block1_pos_modulated_vortex = modulated_vortices{1};
        lh_block2_pos_modulated_vortex = modulated_vortices{2};

        [spontaneuos_vortices, spontaneuous_vortex_longest_duration, modulated_vortices, modulated_vortices_longest_duration] = task_onset_id(spiral_filt_nega,spiral_filt_nega_centreONLY_real, lh_cues, [block_1; block_2]);


        lh_block1_neg_longest_spontaneous = spontaneuous_vortex_longest_duration(1);
        lh_block2_neg_longest_spontaneous = spontaneuous_vortex_longest_duration(2);

        lh_block1_neg_longest_modulated = modulated_vortices_longest_duration(1);
        lh_block2_neg_longest_modulated = modulated_vortices_longest_duration(2);


        lh_block1_nega_spontaneous_vortex = spontaneuos_vortices{1};
        lh_block2_nega_spontaneous_vortex = spontaneuos_vortices{2};

        lh_block1_nega_modulated_vortex = modulated_vortices{1};
        lh_block2_nega_modulated_vortex = modulated_vortices{2};

    end
    
    % first movement block
    block1_rightHand_pos_neg{1,1} = rh_block1_pos_modulated_vortex;
    block1_rightHand_pos_neg{1,2} = rh_block1_pos_longest_modulated;
    block1_rightHand_pos_neg{1,3} = rh_block1_pos_spontaneous_vortex;
    block1_rightHand_pos_neg{1,4} = rh_block1_pos_longest_spontaneous;
    block1_rightHand_pos_neg{1,5} = rh_block1_nega_modulated_vortex;
    block1_rightHand_pos_neg{1,6} = rh_block1_neg_longest_modulated;
    block1_rightHand_pos_neg{1,7} = rh_block1_nega_spontaneous_vortex;
    block1_rightHand_pos_neg{1,8} = rh_block1_neg_longest_spontaneous;

    
    block1_leftHand_pos_neg{1,1} = lh_block1_pos_modulated_vortex;
    block1_leftHand_pos_neg{1,2} = lh_block1_pos_longest_modulated;
    block1_leftHand_pos_neg{1,3} = lh_block1_pos_spontaneous_vortex;
    block1_leftHand_pos_neg{1,4} = lh_block1_pos_longest_spontaneous;
    block1_leftHand_pos_neg{1,5} = lh_block1_nega_modulated_vortex;
    block1_leftHand_pos_neg{1,6} = lh_block1_neg_longest_modulated;
    block1_leftHand_pos_neg{1,7} = lh_block1_nega_spontaneous_vortex;
    block1_leftHand_pos_neg{1,8} = lh_block1_neg_longest_spontaneous;


    block1_rightFoot_pos_neg{1,1} = rf_block1_pos_modulated_vortex;
    block1_rightFoot_pos_neg{1,2} = rf_block1_pos_longest_modulated;
    block1_rightFoot_pos_neg{1,3} = rf_block1_pos_spontaneous_vortex;
    block1_rightFoot_pos_neg{1,4} = rf_block1_pos_longest_spontaneous;
    block1_rightFoot_pos_neg{1,5} = rf_block1_nega_modulated_vortex;
    block1_rightFoot_pos_neg{1,6} = rf_block1_neg_longest_modulated;
    block1_rightFoot_pos_neg{1,7} = rf_block1_nega_spontaneous_vortex;
    block1_rightFoot_pos_neg{1,8} = rf_block1_neg_longest_spontaneous;


    block1_leftFoot_pos_neg{1,1} = lf_block1_pos_modulated_vortex;
    block1_leftFoot_pos_neg{1,2} = lf_block1_pos_longest_modulated;
    block1_leftFoot_pos_neg{1,3} = lf_block1_pos_spontaneous_vortex;
    block1_leftFoot_pos_neg{1,4} = lf_block1_pos_longest_spontaneous;
    block1_leftFoot_pos_neg{1,5} = lf_block1_nega_modulated_vortex;
    block1_leftFoot_pos_neg{1,6} = lf_block1_neg_longest_modulated;
    block1_leftFoot_pos_neg{1,7} = lf_block1_nega_spontaneous_vortex;
    block1_leftFoot_pos_neg{1,8} = lf_block1_neg_longest_spontaneous;


    block1_tongue_pos_neg{1,1} = t_block1_pos_modulated_vortex;
    block1_tongue_pos_neg{1,2} = t_block1_pos_longest_modulated;
    block1_tongue_pos_neg{1,3} = t_block1_pos_spontaneous_vortex;
    block1_tongue_pos_neg{1,4} = t_block1_pos_longest_spontaneous;
    block1_tongue_pos_neg{1,5} = t_block1_nega_modulated_vortex;
    block1_tongue_pos_neg{1,6} = t_block1_neg_longest_modulated;
    block1_tongue_pos_neg{1,7} = t_block1_nega_spontaneous_vortex;
    block1_tongue_pos_neg{1,8} = t_block1_neg_longest_spontaneous;
    
    % Second Movement block

    block2_rightHand_pos_neg{1,1} = rh_block2_pos_modulated_vortex;
    block2_rightHand_pos_neg{1,2} = rh_block2_pos_longest_modulated;
    block2_rightHand_pos_neg{1,3} = rh_block2_pos_spontaneous_vortex;
    block2_rightHand_pos_neg{1,4} = rh_block2_pos_longest_spontaneous;
    block2_rightHand_pos_neg{1,5} = rh_block2_nega_modulated_vortex;
    block2_rightHand_pos_neg{1,6} = rh_block2_neg_longest_modulated;
    block2_rightHand_pos_neg{1,7} = rh_block2_nega_spontaneous_vortex;
    block2_rightHand_pos_neg{1,8} = rh_block2_neg_longest_spontaneous;

    
    block2_leftHand_pos_neg{1,1} = lh_block2_pos_modulated_vortex;
    block2_leftHand_pos_neg{1,2} = lh_block2_pos_longest_modulated;
    block2_leftHand_pos_neg{1,3} = lh_block2_pos_spontaneous_vortex;
    block2_leftHand_pos_neg{1,4} = lh_block2_pos_longest_spontaneous;
    block2_leftHand_pos_neg{1,5} = lh_block2_nega_modulated_vortex;
    block2_leftHand_pos_neg{1,6} = lh_block2_neg_longest_modulated;
    block2_leftHand_pos_neg{1,7} = lh_block2_nega_spontaneous_vortex;
    block2_leftHand_pos_neg{1,8} = lh_block2_neg_longest_spontaneous;


    block2_rightFoot_pos_neg{1,1} = rf_block2_pos_modulated_vortex;
    block2_rightFoot_pos_neg{1,2} = rf_block2_pos_longest_modulated;
    block2_rightFoot_pos_neg{1,3} = rf_block2_pos_spontaneous_vortex;
    block2_rightFoot_pos_neg{1,4} = rf_block2_pos_longest_spontaneous;
    block2_rightFoot_pos_neg{1,5} = rf_block2_nega_modulated_vortex;
    block2_rightFoot_pos_neg{1,6} = rf_block2_neg_longest_modulated;
    block2_rightFoot_pos_neg{1,7} = rf_block2_nega_spontaneous_vortex;
    block2_rightFoot_pos_neg{1,8} = rf_block2_neg_longest_spontaneous;


    block2_leftFoot_pos_neg{1,1} = lf_block2_pos_modulated_vortex;
    block2_leftFoot_pos_neg{1,2} = lf_block2_pos_longest_modulated;
    block2_leftFoot_pos_neg{1,3} = lf_block2_pos_spontaneous_vortex;
    block2_leftFoot_pos_neg{1,4} = lf_block2_pos_longest_spontaneous;
    block2_leftFoot_pos_neg{1,5} = lf_block2_nega_modulated_vortex;
    block2_leftFoot_pos_neg{1,6} = lf_block2_neg_longest_modulated;
    block2_leftFoot_pos_neg{1,7} = lf_block2_nega_spontaneous_vortex;
    block2_leftFoot_pos_neg{1,8} = lf_block2_neg_longest_spontaneous;


    block2_tongue_pos_neg{1,1} = t_block2_pos_modulated_vortex;
    block2_tongue_pos_neg{1,2} = t_block2_pos_longest_modulated;
    block2_tongue_pos_neg{1,3} = t_block2_pos_spontaneous_vortex;
    block2_tongue_pos_neg{1,4} = t_block2_pos_longest_spontaneous;
    block2_tongue_pos_neg{1,5} = t_block2_nega_modulated_vortex;
    block2_tongue_pos_neg{1,6} = t_block2_neg_longest_modulated;
    block2_tongue_pos_neg{1,7} = t_block2_nega_spontaneous_vortex;
    block2_tongue_pos_neg{1,8} = t_block2_neg_longest_spontaneous;        

    if flagSur == 0
        foldername = [main_folder,'/Sample Data/Motor Task/Analysis/'];
        filename = ['task_specific_spiral_ID_leftHem_sub',num2str(subject),'.mat'];
        save([foldername,filename],'block1_rightHand_pos_neg','block1_leftHand_pos_neg','block1_rightFoot_pos_neg','block1_leftFoot_pos_neg','block1_tongue_pos_neg',...
                                   'block2_rightHand_pos_neg','block2_leftHand_pos_neg','block2_rightFoot_pos_neg','block2_leftFoot_pos_neg','block2_tongue_pos_neg');
    elseif flagSur == 1
        % foldername = [main_folder,'/Sample Data/Motor Task/Analysis/'];
        % filename = ['task_specific_spiral_distribution_leftnew_sur.mat'];
        % save([foldername,filename],'RH_block1_pos_nega_VOI','LH_block1_pos_nega_VOI','RF_block1_pos_nega_VOI','LF_block1_pos_nega_VOI','T_block1_pos_nega_VOI','RH_block2_pos_nega_VOI','LH_block2_pos_nega_VOI','RF_block2_pos_nega_VOI','LF_block2_pos_nega_VOI','T_block2_pos_nega_VOI');

    end
    cd(main_folder)

end

end
