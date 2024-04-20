function [spontaneuos_vortices, spontaneuous_vortex_longest_duration, modulated_vortices, modulated_vortices_longest_duration] = task_onset_id(spiral_filt,centres, cues, blocks)

block1_spontaneous_longest_duration = 0;
block1_modulated_longest_duration = 0;
block2_spontaneous_longest_duration = 0;
vortexIndex_block1 = 0;
vortexIndex_block2 = 0;
modulatedIndex_block1 = 0;
modulatedIndex_block2 = 0;
possible_modulation = 0;

for ipatt = 1:size(spiral_filt,1)
    for time = blocks(1,:)
        if isempty(spiral_filt{ipatt,time-1}) && ~isempty(spiral_filt{ipatt,time})
            possible_VOI = ipatt;
            VOI_time = find(~cellfun('isempty', spiral_filt(possible_VOI,:)));
            if size(VOI_time,2) > 3 && ismember(VOI_time(1),blocks(1,:))
                vortexIndex_block1 = vortexIndex_block1 + 1;
                block1_vortex_of_interest{vortexIndex_block1,1} = possible_VOI;
                block1_vortex_of_interest{vortexIndex_block1,2} = centres(ipatt,VOI_time);
                block1_vortex_of_interest{vortexIndex_block1,3} = VOI_time;

                if VOI_time(end) > block1_spontaneous_longest_duration
                    block1_spontaneous_longest_duration = VOI_time(end);
                end
            end
        elseif ~isempty(spiral_filt{ipatt,blocks(1,1)-1}) && ~isempty(spiral_filt{ipatt,blocks(1,1)}) && ipatt ~= possible_modulation
            possible_modulation = ipatt;
            modulation_time = find(~cellfun('isempty', spiral_filt(possible_modulation,:)));
            if size(modulation_time,2) > 3 && modulation_time(1) < blocks(1,1)     %cues(1)   %Cue time gonna consider only less than cue (spontaneuos spirals between cue and onset may be something else)
                modulatedIndex_block1 = modulatedIndex_block1 + 1;
                block1_modulated_vortices{modulatedIndex_block1,1} = possible_modulation;
                block1_modulated_vortices{modulatedIndex_block1,2} = centres(ipatt,modulation_time);
                block1_modulated_vortices{modulatedIndex_block1,3} = modulation_time;
                if modulation_time(end) > block1_modulated_longest_duration
                    block1_modulated_longest_duration = modulation_time(end);
                end
            end
        end
    end

    for time = blocks(2,:)
        if isempty(spiral_filt{ipatt,time-1}) && ~isempty(spiral_filt{ipatt,time})
            possible_VOI = ipatt;
            VOI_time = find(~cellfun('isempty', spiral_filt(possible_VOI,:)));
            if size(VOI_time,2) > 3 && ismember(VOI_time(1),blocks(2,:))
                vortexIndex_block2 = vortexIndex_block2 + 1;
                block2_vortex_of_interest{vortexIndex_block2,1} = possible_VOI;
                block2_vortex_of_interest{vortexIndex_block2,2} = centres(ipatt,VOI_time);
                block2_vortex_of_interest{vortexIndex_block2,3} = VOI_time;
                if VOI_time(end) >block2_spontaneous_longest_duration
                    block2_spontaneous_longest_duration = VOI_time(end);
                end
            end
        elseif ~isempty(spiral_filt{ipatt,blocks(2,1)-1}) && ~isempty(spiral_filt{ipatt,blocks(2,1)}) && ipatt ~= possible_modulation
            possible_modulation = ipatt;
            modulation_time = find(~cellfun('isempty', spiral_filt(possible_modulation,:)));
            if size(modulation_time,2) > 3 && modulation_time(1) < blocks(2,1)   %Cue time gonna consider only less than cue (spontaneuos spirals between cue and onset may be something else)
                modulatedIndex_block2 = modulatedIndex_block2 +1;
                block2_modulated_vortices{modulatedIndex_block2,1} = possible_modulation;
                block2_modulated_vortices{modulatedIndex_block2,2} = centres(ipatt,modulation_time);
                block2_modulated_vortices{modulatedIndex_block2,3} = modulation_time;
                if modulation_time(end) > block1_modulated_longest_duration
                    block2_modulated_longest_duration = modulation_time(end);
                end
            end
        end
    end
end

spontaneuos_vortices{1} = block1_vortex_of_interest;
spontaneuos_vortices{2} = block2_vortex_of_interest;
modulated_vortices{1} = block1_modulated_vortices ;
modulated_vortices{2} = block2_modulated_vortices;

spontaneuous_vortex_longest_duration = [block1_spontaneous_longest_duration block2_spontaneous_longest_duration];
modulated_vortices_longest_duration = [block1_modulated_longest_duration block2_modulated_longest_duration];


end