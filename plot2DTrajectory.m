
    folders = get_folders(rootDir);
   
    %% Process each folder
    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
for ff = 1:folderNum
      
      %% Get folder information
  folder = folders(ff).folder;
      
       
%% Load in fictrac & ROI data
    if strcmp(folder(end),'.')
        folder = folder(1:end-2); 
    end

    % Get data files
    expID = get_expID(folder);
    expList = {expID};
    
    % Load metadata 
    [~, trialMd] = load_metadata(expList, folder);

    % Load imaging data
    roiData = load_roi_data(expList, folder);
    
    % Load FicTrac data
    [~,~,ftData_dat] = load_ft_data(expList, folder, 0, 1);



    numTrials = max(size(unique(roiData.trialNum),1),length(trialMd.trialNum)); 
    %%
    for nTrial = 1:numTrials

        figure();patch(ftData_dat.intX{1},ftData_dat.intY{1},ftData_dat.trialTime{1},'EdgeColor','interp','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
        a=colorbar; a.Label.String = 'Time (sec)';
        set(gcf,'color','w');
        axis off

        saveas(gcf, fullfile(folder,'images',[expID,'_',num2str(nTrial),'_pathTrajectory.fig']));
        saveas(gcf, fullfile(folder,'images',[expID,'_',num2str(nTrial),'_pathTrajectory.png']));
        %saveas(gcf, fullfile(folder,'images',[expID,'_',num2str(nTrial),'_pathTrajectory.svg']));
    end
end

close all 