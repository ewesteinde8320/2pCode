function plot_brainSpace_activity_byHeading(folder, goal,rho) 

window = 90; 

threshold=3;
trial_count = 1;

trial_vy = [];
trial_angle = [];
trial_vf = []; 
speed_sum = [];    

expID = get_expID(folder);
expList = {expID};

nTrial = 1; 

[~,ftT, ~] = load_ft_data(expList, folder, 1, 0);

% Load metadata 
[expMd, trialMd] = load_metadata(expList, folder);

% Load imaging data
roiData = load_roi_data(expList, folder);

processedData_dir = fullfile(folder,'processed_data');

data_filelist = dir(processedData_dir);
for files = 1:length(data_filelist)
    if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
        load(fullfile(processedData_dir,data_filelist(files).name));
    end
    load(fullfile(processedData_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']))
end

total_mov_mm = abs(ftT.velFor{1}) + abs(ftT.velSide{1}) + abs(ftT.velYaw{1}*4.5);
no0vel_idx = find(total_mov_mm > threshold);
angle = -ftT.cueAngle{1};
angle = angle(no0vel_idx); 
angle = wrapTo180(wrapTo180(angle)-wrapTo180((goal)));
edges_angle = [-180 - window/2:window:180 + window/2]; 

Z = []; 
for roi = 1:size(ZData,1)
    Z(:,roi) = ZData.(3){roi}(no0vel_idx); 
end

roiActivity = []; 
    for roi = 1:size(Z,2)
        activity = Z(:,roi);
        if isnan(activity) | isnan(angle)
            activity = activity( ~isnan(activity) & ~isnan(angle)); 
            angle = angle( ~isnan(activity) & ~isnan(angle)); 
        end
        [zscore, centers_angle] = binData(activity, angle, edges_angle);
        roiActivity(roi,:) = zscore; 
    end
    
    roiActivity(:,1) = (roiActivity(:,1) + roiActivity(:,end))./2;
    roiActivity(:,end) = []; 
    
    %[bump_params, ~, ~] = fit_sinusoid(roiActivity,[0,2*pi], 0);

%%
    activityTable = Z;
    ncol = size(roiActivity,2);
    color = (cbrewer2('Greens', ncol));
    
    test = circshift(color,8,1);
    
    brainSpace = linspace(180,-180,size(Z,2));
    legendNames = []; 
    legendNames(:,1)= centers_angle(1:end-1);
    legendNames = num2str(legendNames);
    
    figure(1);clf;
    set(gcf,'color','w','renderer','painters')
    hold on
    for h = 1:size(roiActivity,2)
        plot(brainSpace,roiActivity(:,h),'color',test(h ,:),'LineWidth',1.5)
    end
    legend({legendNames})
    xlabel('brain space')
    ylabel('z-scored df/f')
    legend boxoff
    title(num2str(rho))


end




