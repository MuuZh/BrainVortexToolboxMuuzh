clear

% read data
filepaths = get_spiral_detected_path("LSD", "LSD");

usingvars = {'significant_spiral_duration_extend';
            'spiral_filt_nega_centreONLY_real';
            'spiral_filt_nega_real';
            'spiral_filt_pos_centreONLY_real';
            'spiral_filt_pos_real';
            'spiral_size_real_95perc_accu'
            };
for i = 1:1
    spiraldetecetd = load(filepaths{i},usingvars{:});
end