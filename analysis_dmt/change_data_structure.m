clear

func = "DMT";
group = "PCB";
% read data
filepaths = get_spiral_detected_path(func, group);

usingvars = {'significant_spiral_duration_extend';
            'spiral_filt_nega_centreONLY_real';
            'spiral_filt_nega_real';
            'spiral_filt_pos_centreONLY_real';
            'spiral_filt_pos_real';
            'spiral_size_real_95perc_accu'
            };


            
for i = 1:length(filepaths)
    spiraldetecetd = load(filepaths{i},usingvars{:});
    
    spiral_pos_matrix_index = sparse2list(spiraldetecetd.spiral_filt_pos_real, 1);
    spiral_neg_matrix_index = sparse2list(spiraldetecetd.spiral_filt_nega_real, -1);

    spiral_pos_center_coord = spiraldetecetd.spiral_filt_pos_centreONLY_real;
    spiral_neg_center_coord = spiraldetecetd.spiral_filt_nega_centreONLY_real;

    significant_spiral_duration_extend = spiraldetecetd.significant_spiral_duration_extend;
    spiral_size_real_95perc_accu = spiraldetecetd.spiral_size_real_95perc_accu;
    
    fn_to_save = sprintf('%s_%s_S%02d.mat', func, group, i);

    SUB_ID = i;

    save(fn_to_save,...
        'spiral_pos_matrix_index', 'spiral_neg_matrix_index', ...
        'spiral_pos_center_coord', 'spiral_neg_center_coord', ...
        'significant_spiral_duration_extend', 'spiral_size_real_95perc_accu', ...
        'func', 'group', 'SUB_ID' ...
        );

end
