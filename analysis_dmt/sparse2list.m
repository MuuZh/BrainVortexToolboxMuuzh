function all_matrix_index = sparse2list(spin_pos_span, spin)
    if spin > 0
        disp('Spin is positive.');
        currentspin = 1;
    elseif spin < 0
        disp('Spin is negative.');
        currentspin = -1;
    elseif spin == 0
        disp('Now just for centers.');
    end
    [spiral_num, frame_num] = size(spin_pos_span);
    all_matrix_index = cell(spiral_num, frame_num);
    
    for spiral_idx = 1:spiral_num
        for frame_idx = 1:frame_num
            if ~isempty(spin_pos_span{spiral_idx, frame_idx})
                currentcell = spin_pos_span{spiral_idx, frame_idx};
                if spin == 0
                    all_matrix_index{spiral_idx, frame_idx} = currentcell;
                else
                    [row, col, val] = find(currentcell);
                    if currentspin == 1
                        indices = val > 0;
                    elseif currentspin == -1
                        indices = val < 0;
                    end
                    selected_rows = row(indices);
                    selected_cols = col(indices);
                    all_matrix_index{spiral_idx, frame_idx} = [selected_rows, selected_cols];
                end
            end
        end
    end

end