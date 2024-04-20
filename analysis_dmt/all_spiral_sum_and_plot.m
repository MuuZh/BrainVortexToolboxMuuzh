clear
close all
orig_folder = pwd;
cd('E:\resarch_data\fMRI\BrainVortexToolbox-main')

func = 'DMT';
% func = 'PCB';

load("parcellation_template.mat")
load("parcellation_template22_RightBrain_subject1-100.mat")


if func == 'DMT'
    subnums = 22;
    fn = ['E:\resarch_data\fMRI\fMRI_DMT\BrainToolBoxed\DMT_post1\forplotting\','alltime_spiral_grid_DMT.mat']
    load(fn);
elseif func == 'PCB'
    subnums = 20;
    fn = ['E:\resarch_data\fMRI\fMRI_DMT\BrainToolBoxed\PCB_post1\forplotting\','alltime_spiral_grid_PCB.mat']
    load(fn);
end

all_spirals_neg;
all_spirals_pos;
all_spirals_overall;



% subplot for 14 subjects
figure('units','normalized','outerposition',[0 0 1 1])

cop = 'jet'
colormap(cop)



for i=1:subnums
    subplot(5,5,i)
    hold on
    set(gca,'ydir','normal')
    imagesc((abs(all_spirals_neg{i}) + abs(all_spirals_pos{i})))
    for parcellation_ID = 1:22 % Numbers correspond to certain areas of brain
        parcellation_template_1par = parcellation_template;
        parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
        
        parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
        
        B = bwboundaries(parcellation_template_1par(:,:,1),'noholes');
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])
            
        end
    end
    axis equal
    xlim([0 251])
    ylim([0 175])
    colorbar()
%     title([func ' Subject ' sprintf('%d',i) ' counting both pos and neg'])
    

end

figure('units','normalized','outerposition',[0 0 1 1])

for i=1:subnums
    if i == 1
        all_spirals_combined = (abs(all_spirals_neg{i}) + abs(all_spirals_pos{i}));
    else
        all_spirals_combined = all_spirals_combined + (abs(all_spirals_neg{i}) + abs(all_spirals_pos{i}));
    end
end

hold on
% set(gca,'ydir','normal')


imagesc((all_spirals_combined))

for parcellation_ID = 1:22 % Numbers correspond to certain areas of brain
    parcellation_template_1par = parcellation_template;
    parcellation_template_1par(isnan(parcellation_template_1par)) = 0;
    
    parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;
    
    B = bwboundaries(parcellation_template_1par(:,:,1),'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])
        
    end
end

axis equal
xlim([0 251])
ylim([0 175])
colormap(cop)
colorbar()
% title([func ' All subjects spirals combined counting both pos and neg'])
ax=gca;
ax.XAxis.Visible='off';
ax.YAxis.Visible='off';
colorbar('FontSize',26);

mean(reshape(all_spirals_combined,[],1))
std(reshape(all_spirals_combined,[],1))
figure
hist(reshape(all_spirals_combined,[],1))

% if func == 'LSD'
%     save('alltime_spiral_grid_LSD_forjupter.mat', 'all_spirals_combined');
% elseif func == 'PCB'
%     save('alltime_spiral_grid_PCB_forjupter.mat', 'all_spirals_combined');
% end

cd(orig_folder)