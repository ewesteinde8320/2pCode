function MakeJumpPFL2Summary_allTrial(rootDir)       

     folders = get_folders(rootDir,1,0);
     
     headers = {'folder','trial','jumpSize','jumpIdx','velFor','velSide','velYaw','cueAngle','corrected','timeMov','Z','bumpAmp','region'};
     jump_summary = cell2table(cell(1,13),'VariableNames',headers); 
     count = 1; 

    %% Process each folder
    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
    countFail = 1;
    for ff = 1:folderNum
      
      %% Get folder information
      folder = folders(ff).folder;
      
       
%% Load in fictrac & ROI data
       try
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end

             processedData_dir = fullfile(folder,'processed_data');

            % Get data files
            expID = get_expID(folder);
            expList = {expID};

            % Load metadata 
            [expMd, trialMd] = load_metadata(expList, folder);

            % Load imaging data
            roiData = load_roi_data(expList, folder);

            % Load FicTrac data
            [~,ftData, ~] = load_ft_data(expList, folder, 1, 0);

            % Load panels metadata
            panelsMetadata = load_panels_metadata(expList, folder);

            try
                numTrials = max(size(unique(roiData.trialNum),1),length(trialMd.trialNum)); 
            catch
                numTrials = 1; 
            end
        %%
        for nTrial = 1:numTrials
        %% get processed data
            bump_params = [];
            data_filelist = dir(processedData_dir);
            for files = 1:length(data_filelist)
                if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                    load(fullfile(processedData_dir,data_filelist(files).name));
                end
            end
            
            load(fullfile(processedData_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']))

    %% Remove idx where the fly isn't moving
            %ftT = ftT_minSmooth; 
            window = 15; 
            [~, jump_array, ~] = detect_jumps(ftT, window, 5,0);
            [~, jump_array_rho, ~] = detect_jumps(ftT, 30, 5,0);

                count90 = 1;
                count180 = 1;
                jump90_time = []; 
                jump180_time = []; 
                all_roi = []; 
                
                    jumpIdx = [jump_array(1,1):jump_array(1,3)]';  
                    jumpIdx_rho = [jump_array_rho(1,1):jump_array_rho(1,3)]';
                    
                    if jumpIdx(1) < 1 || jumpIdx_rho(1) < 1
                        jump_array(1,:) = []; 
                        jump_array_rho(1,:) = []; 
                    end

                    for jump = 1:size(jump_array,1)
                        jumpIdx = [jump_array(jump,1):jump_array(jump,3)]'; 
                        if jumpIdx(end) > size(ftT.velFor{1},1)
                            jump_array(jump,:) = []; 
                            jump_array_rho(jump,:) = []; 
                        end
                    end 
                    
                for jump = 1:size(jump_array,1)
                    jumpSize = jump_array(jump,4);
                    jumpIdx = [jump_array(jump,1):jump_array(jump,3)]';  
                    jump_vf = ftT.velFor{1}(jumpIdx);
                    jump_vy = (rad2deg(ftT.velYaw{1}(jumpIdx)));
                    jump_vs = (ftT.velSide{1}(jumpIdx));
                    jump_angle = ftT.cueAngle{1}(jumpIdx); 
                    
                    for roi = 1:size(ZData,1)
                        all_roi(roi,:) = ZData.Z{roi}(jumpIdx);
                    end
                    
                    jump_Z = mean(all_roi,1,'omitnan');
                    jump_bumpAmp = bump_params.amp(jumpIdx);
                    
                    preJumpIdx = [jump_array(jump,1):jump_array(jump,2)-1];
                    preJumpIdx_rho = [jump_array_rho(jump,1):jump_array_rho(jump,2)-1];
                    timeMov = (sum(abs(ftT.velFor{1}(preJumpIdx_rho)) + abs(ftT.velSide{1}(preJumpIdx_rho)) + abs(deg2rad(ftT.velYaw{1}(preJumpIdx_rho))*4.5) > 1.5))/60; % seconds
                    [rho, ~] = CalculateAverageHeading(ftT,1.5, preJumpIdx_rho);
                    [~, preJumpAngle] = CalculateAverageHeading(ftT,0, preJumpIdx);
                    
                    %preJumpAngle = ftT.cueAngle{1}(preJumpIdx(end));
                    
                    postJumpIdx = jump_array(jump,2)+1;
                    sizePostJump = size(ftT.cueAngle{1}(postJumpIdx:jump_array(jump,3)));
                    preJump = ones(sizePostJump) * preJumpAngle;
                    cueDiff = angdiff(preJump,deg2rad(ftT.cueAngle{1}(postJumpIdx:jump_array(jump,3)))); 
                    
                    if jumpSize == 180
                        correct = sum(abs(cueDiff) < deg2rad(45)); 
                    else
                        correct = sum(abs(cueDiff) < deg2rad(30));
                    end
                    
                    jump_summary.folder(count) = {folder}; 
                    jump_summary.trial(count) = {nTrial}; 
                    jump_summary.jumpSize(count) = {jumpSize}; 
                    jump_summary.jumpIdx(count) = {jumpIdx};
                    jump_summary.velFor(count) = {jump_vf};
                    jump_summary.velYaw(count) = {jump_vy};
                    jump_summary.velSide(count) = {jump_vs}; 
                    jump_summary.cueAngle(count) = {jump_angle};
                    jump_summary.timeMov(count) = {timeMov}; 
                    jump_summary.Z(count) = {jump_Z}; 
                    jump_summary.bumpAmp(count) = {jump_bumpAmp}; 
                    isPB = regexp(folder,'PB');
                    isFB = regexp(folder,'FB');
                    
                    if ~isempty(isPB)
                        jump_summary.region(count) = {'PB'};
                    elseif ~isempty(isFB)
                        jump_summary.region(count) = {'FB'};
                    else
                        jump_summary.region(count) = {'LAL'};
                    end

                  
                    if  timeMov >= 5  && correct >= 1 %&& rho >= 0.88
                        jump_summary.corrected(count) = {1};
             
%                     elseif correct < 1 && timeMov >= 5 && rho >= 0.88
%                         jump_summary.corrected(count) = {2}; 
%                         figure(1);clf;plot(jump_angle)
                    elseif correct < 1 && timeMov >= 5
                        jump_summary.corrected(count) = {0}; 

                    else
                        jump_summary.corrected(count) = {3}; 
                    
                        % corrected, rho > threshold, timeMov < threshold
                    end
                    
                    count = count + 1; 
                    
                end
        end
        catch
            disp([folder, ' failed'])
        end
    end  
    %%
    Sum180 = jump_summary(cell2mat(jump_summary.jumpSize) == 180 & cell2mat(jump_summary.corrected) == 1,:); 
    Sum180_not = jump_summary(cell2mat(jump_summary.jumpSize) == 180 & cell2mat(jump_summary.corrected) == 0,:); 
    Sum90 = jump_summary((cell2mat(jump_summary.jumpSize)) == 90& cell2mat(jump_summary.corrected) == 1,:);
    Summ90 = jump_summary(cell2mat(jump_summary.jumpSize) == -90 & cell2mat(jump_summary.corrected) == 0,:);
%%
    jumpTime = [-window+(1/60):1/60:window+(1/60)]; 
    ZSum = zeros(size(Sum180,1),size(jump_vf,1));
    bAmpSum = zeros(size(Sum180,1),size(jump_vf,1));
    vfSum = zeros(size(Sum180,1),size(jump_vf,1));
    vySum = zeros(size(Sum180,1),size(jump_vf,1));
    vsSum = zeros(size(Sum180,1),size(jump_vf,1));
    figure();
    sgtitle('180 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    ax1 = subplot(1,1,1); 
    hold on
%     ax2 = subplot(4,1,3);
%     hold on
%     ax3 = subplot(4,1,4); 
%     hold on
%     ax4 = subplot(4,1,2); 
%     hold on

    for jump = 1:size(Sum180,1)  
%         figure();
%         h(1) = subplot(3,1,1);
%         plot( jumpTime,smoothdata(abs(Sum180.velYaw{jump})','loess',5))
%         h(2) = subplot(3,1,2);
%         plot( jumpTime,smoothdata(Sum180.bumpAmp{jump},'loess',5))
%         h(3) = subplot(3,1,3);
%         plot(jumpTime, Sum180.cueAngle{jump})
%         linkaxes(h,'x')
        ZSum(jump,:) = (Sum180.Z{jump});
        bAmpSum(jump,:) = (Sum180.bumpAmp{jump});
        vfSum(jump,:) = smoothdata(Sum180.velFor{jump}','loess',50); 
        vySum(jump,:) = abs(smoothdata(Sum180.velYaw{jump}','loess',50));
        vsSum(jump,:) = smoothdata(Sum180.velSide{jump}','loess',50);
    end
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_bAmp = std(bAmpSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum(:,keepIndex),'omitnan') + SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum(:,keepIndex),'omitnan') - SEM_Z(keepIndex)];
    plot(ax1,jumpTime,mean(ZSum,'omitnan'),'k','LineWidth',1.5)
    patch(ax1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(ax1,mean(smoothdata(mean(ZSum,'omitnan'),'gaussian',75)),'k')
    ylabel(ax1,'Z')
    xlabel(ax1,'Time (s)')
    xlim(ax1,[min(jumpTime),max(jumpTime)])
    ylim(ax1,[-0.2,0.8])
    
%     keepIndex = ~isnan(SEM_bAmp);
%     SEMhigh = [mean(bAmpSum(:,keepIndex),'omitnan') + SEM_bAmp(keepIndex)]; 
%     SEMlow = [mean(bAmpSum(:,keepIndex),'omitnan') - SEM_bAmp(keepIndex)];
%     plot(ax4,jumpTime,smoothdata(mean(bAmpSum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
%     patch(ax4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     %yline(cx1,mean(smoothdata(mean(RLALSum,'omitnan'),'gaussian',75)),'k')
%     ylim(ax4,[0.6,1.8])
%     xlim(ax4,[min(jumpTime),max(jumpTime)])
%     ylabel(ax4,'bump amp')
%     xlabel(ax4,'Time (s)')
%     
%     
%     keepIndex = ~isnan(SEM_vf);
%     SEMhigh = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vf(keepIndex)]; 
%     SEMlow = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vf(keepIndex)];
%     plot(ax2,jumpTime,smoothdata(mean(vfSum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
%     patch(ax2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     %yline(ax2,mean(smoothdata(mean(vfSum,'omitnan'),'gaussian',75)),'k')
%     ylim(ax2,[0,4])
%     xlim(ax2,[min(jumpTime),max(jumpTime)])
%     ylabel(ax2,'forward velocity (mm/s)')
%     xlabel(ax2,'Time (s)')
%     
%     keepIndex = ~isnan(SEM_vy);
%     SEMhigh = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vy(keepIndex)]; 
%     SEMlow = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vy(keepIndex)];
%     plot(ax3,jumpTime,smoothdata(mean(vySum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
%     patch(ax3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     %yline(ax3,mean(smoothdata(mean(vySum,'omitnan'),'gaussian',75)),'k')
%     ylim(ax3,[15,80])
%     xlim(ax3,[min(jumpTime),max(jumpTime)])
%     ylabel(ax3,'rotational speed (deg/s)')
%     xlabel(ax3,'Time (s)')
    
%     %
    jumpTime = [-window+(1/60):1/60:window+(1/60)]; 
    ZSum = zeros(size(Sum180_not,1),size(jump_vf,1));
    bAmpSum = zeros(size(Sum180_not,1),size(jump_vf,1));
    vfSum = zeros(size(Sum180_not,1),size(jump_vf,1));
    vySum = zeros(size(Sum180_not,1),size(jump_vf,1));
    vsSum = zeros(size(Sum180_not,1),size(jump_vf,1));

    for jump = 1:size(Sum180_not,1)  
        ZSum(jump,:) = smoothdata(Sum180_not.Z{jump},'loess',100);
        bAmpSum(jump,:) = smoothdata(Sum180_not.bumpAmp{jump},'loess',100);
        vfSum(jump,:) = smoothdata(Sum180_not.velFor{jump}','loess',100); 
        vySum(jump,:) = smoothdata(abs(Sum180_not.velYaw{jump})','loess',100);
        vsSum(jump,:) = abs(Sum180_not.velSide{jump})';
    end
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_bAmp = std(bAmpSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum(:,keepIndex),'omitnan') + SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum(:,keepIndex),'omitnan') - SEM_Z(keepIndex)];
    plot(ax1,jumpTime,mean(ZSum,'omitnan'),'b','LineWidth',1.5)
    patch(ax1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(ax1,mean(smoothdata(mean(ZSum,'omitnan'),'gaussian',75)),'k')
    ylabel(ax1,'Z')
    xlabel(ax1,'Time (s)')
    xlim(ax1,[min(jumpTime),max(jumpTime)])
%     
%     keepIndex = ~isnan(SEM_bAmp);
%     SEMhigh = [mean(bAmpSum(:,keepIndex),'omitnan') + SEM_bAmp(keepIndex)]; 
%     SEMlow = [mean(bAmpSum(:,keepIndex),'omitnan') - SEM_bAmp(keepIndex)];
%     plot(ax4,jumpTime,smoothdata(mean(bAmpSum,'omitnan'),'gaussian',5),'b','LineWidth',1.5)
%     patch(ax4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     %yline(cx1,mean(smoothdata(mean(RLALSum,'omitnan'),'gaussian',75)),'k')
%     %ylim(cx1,[-0.92,0.55])
%     xlim(ax4,[min(jumpTime),max(jumpTime)])
%     ylabel(ax4,'bump amp')
%     xlabel(ax4,'Time (s)')
%     
%     
%     keepIndex = ~isnan(SEM_vf);
%     SEMhigh = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vf(keepIndex)]; 
%     SEMlow = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vf(keepIndex)];
%     plot(ax2,jumpTime,smoothdata(mean(vfSum,'omitnan'),'gaussian',5),'b','LineWidth',1.5)
%     patch(ax2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0],'FaceAlpha',0.25,'EdgeColor','none')
%     %yline(ax2,mean(smoothdata(mean(vfSum,'omitnan'),'gaussian',75)),'k')
%     ylim(ax2,[0,3.5])
%     xlim(ax2,[min(jumpTime),max(jumpTime)])
%     ylabel(ax2,'forward velocity (mm/s)')
%     xlabel(ax2,'Time (s)')
%     
%     keepIndex = ~isnan(SEM_vy);
%     SEMhigh = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vy(keepIndex)]; 
%     SEMlow = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vy(keepIndex)];
%     plot(ax3,jumpTime,smoothdata(mean(vySum,'omitnan'),'gaussian',5),'b','LineWidth',1.5)
%     patch(ax3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     %yline(ax3,mean(smoothdata(mean(vySum,'omitnan'),'gaussian',75)),'k')
%     ylim(ax3,[0,100])
%     xlim(ax3,[min(jumpTime),max(jumpTime)])
%     ylabel(ax3,'rotational speed (deg/s)')
%     xlabel(ax3,'Time (s)')
    
    %

     linkaxes([ax1,ax2,ax3,ax4],'x')
        
    vfSum = zeros(size(Sum90,1),size(jump_vf,1));
    vySum = zeros(size(Sum90,1),size(jump_vf,1));
    vsSum = zeros(size(Sum90,1),size(jump_vf,1));
    ZSum = zeros(size(Sum90,1),size(jump_vf,1));
    bAmpSum = zeros(size(Sum90,1),size(jump_vf,1));
    figure();
    sgtitle('90 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    bx1 = subplot(4,1,1);
    hold on
    bx2 = subplot(4,1,3);
    hold on
    bx3 = subplot(4,1,4); 
    hold on
    bx4 = subplot(4,1,2); 
    hold on
    for jump = 1:size(Sum90,1)  
        ZSum(jump,:) = Sum90.Z{jump};
        bAmpSum(jump,:) = Sum90.bumpAmp{jump};
        vfSum(jump,:) = smoothdata(Sum90.velFor{jump}','loess',50); 
        vySum(jump,:) = smoothdata(abs(Sum90.velYaw{jump})','loess',50);
        vsSum(jump,:) = smoothdata(Sum90.velSide{jump}','loess',50);
    end
    
        startCount = jump;
    for jump = 1:size(Summ90,1)
        ZSum(jump,:) = Summ90.Z{jump};
        bAmpSum(jump,:) = Summ90.bumpAmp{jump};
        vfSum(startCount,:) = smoothdata(Summ90.velFor{jump}','loess',100); 
        %plot(bx3,jumpTime,Sum90.velYaw{jump})
        vySum(startCount,:) = smoothdata(abs(Summ90.velYaw{jump})','loess',100);
        %plot(bx4,jumpTime,Sum90.velSide{jump})
        vsSum(startCount,:) = smoothdata(Summ90.velSide{jump}','loess',100);
        startCount = startCount + 1; 
    end
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vs = std(vsSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_bAmp = std(bAmpSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum(:,keepIndex),'omitnan') + SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum(:,keepIndex),'omitnan') - SEM_Z(keepIndex)];
    plot(bx1,jumpTime,mean(ZSum,'omitnan'),'k','LineWidth',1.5)
    patch(bx1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(bx1,0,'k')
    ylabel(bx1,'Z')
    xlabel(bx1,'Time (s)')
    xlim(bx1,[min(jumpTime),max(jumpTime)])
    ylim(bx1,[-0.2,0.8])
    
    keepIndex = ~isnan(SEM_bAmp);
    SEMhigh = [mean(bAmpSum(:,keepIndex),'omitnan') + SEM_bAmp(keepIndex)]; 
    SEMlow = [mean(bAmpSum(:,keepIndex),'omitnan') - SEM_bAmp(keepIndex)];
    plot(bx4,jumpTime,smoothdata(mean(bAmpSum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
    patch(bx4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(cx1,mean(smoothdata(mean(RLALSum,'omitnan'),'gaussian',75)),'k')
    ylim(bx4,[0.6,1.8])
    xlim(bx4,[min(jumpTime),max(jumpTime)])
    ylabel(bx4,'bump amp')
    xlabel(bx4,'Time (s)')
    
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [smoothdata(mean(vfSum(:,keepIndex)),'gaussian',5) + SEM_vf(keepIndex)]; 
    SEMlow = [smoothdata(mean(vfSum(:,keepIndex)),'gaussian',5) - SEM_vf(keepIndex)];
    %xlim(bx1,[min(jumpTime),max(jumpTime)])
    plot(bx2,jumpTime,smoothdata(mean(vfSum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
    patch(bx2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(bx2,mean(smoothdata(mean(vfSum,'omitnan'),'gaussian',5)),'k')
    ylim(bx2,[0,4])
    xlim(bx2,[min(jumpTime),max(jumpTime)])
    ylabel(bx2,'forward velocity (mm/s)')
    xlabel(bx2,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vy(keepIndex)]; 
    SEMlow = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vy(keepIndex)];
    plot(bx3,jumpTime,smoothdata(mean(vySum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
    patch(bx3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(bx3,0,'k')
    ylim(bx3,[[20,80]])
    xlim(bx3,[min(jumpTime),max(jumpTime)])
    ylabel(bx3,'rotational velocity (deg/s)')
    xlabel(bx3,'Time (s)')
    
%     keepIndex = ~isnan(SEM_vs);
%     SEMhigh = [smoothdata(mean(vsSum(:,keepIndex),'omitnan'),'gaussian',75) + SEM_vs(keepIndex)]; 
%     SEMlow = [smoothdata(mean(vsSum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vs(keepIndex)];
%     plot(bx4,jumpTime,smoothdata(mean(vsSum,'omitnan'),'gaussian',75),'k','LineWidth',1.5)
%     patch(bx4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     yline(bx4,0,'k')
%     ylim(bx4,[-1,1])
%     xlim(bx4,[min(jumpTime),max(jumpTime)])
     linkaxes([bx1,bx2,bx3,bx4],'x')
%    
    vfSum = zeros(size(Summ90,1),size(jump_vf,1));
    vySum = zeros(size(Summ90,1),size(jump_vf,1));
    vsSum = zeros(size(Summ90,1),size(jump_vf,1)); 
    ZSum = zeros(size(Summ90,1),size(jump_vf,1));
    bAmpSum = zeros(size(Summ90,1),size(jump_vf,1));
    %RLALSum = zeros(size(Summ90,1),size(jump_vf,1));
    figure();
    sgtitle('-90 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    cx1 = subplot(4,1,1);
    hold on
    cx2 = subplot(4,1,3);
    hold on
    cx3 = subplot(4,1,4); 
    hold on
    cx4 = subplot(4,1,2); 
    hold on
    for jump = 1:size(Summ90,1)  
        ZSum(jump,:) = Summ90.Z{jump};
        bAmpSum(jump,:) = Summ90.bumpAmp{jump};
        vfSum(jump,:) = smoothdata(Summ90.velFor{jump}','loess',50); 
        vySum(jump,:) = smoothdata(Summ90.velYaw{jump}','loess',50);
        vsSum(jump,:) = smoothdata(Summ90.velSide{jump}','loess',50);
    end
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vs = std(vsSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_Z = std(ZSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_bAmp = std(bAmpSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    keepIndex = ~isnan(SEM_Z);
    SEMhigh = [mean(ZSum(:,keepIndex),'omitnan')+ SEM_Z(keepIndex)]; 
    SEMlow = [mean(ZSum(:,keepIndex),'omitnan') - SEM_Z(keepIndex)];
    plot(cx1,jumpTime,mean(ZSum,'omitnan'),'k','LineWidth',1.5)
    patch(cx1,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(cx1,0,'k')
    ylabel(cx1,'Z')
    xlabel(cx1,'Time (s)')
    xlim(cx1,[min(jumpTime),max(jumpTime)])
    ylim(cx1,[-0.2,0.8])
    
    keepIndex = ~isnan(SEM_bAmp);
    SEMhigh = [mean(bAmpSum(:,keepIndex),'omitnan') + SEM_bAmp(keepIndex)]; 
    SEMlow = [mean(bAmpSum(:,keepIndex),'omitnan') - SEM_bAmp(keepIndex)];
    plot(cx4,jumpTime,smoothdata(mean(bAmpSum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
    patch(cx4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(cx1,mean(smoothdata(mean(RLALSum,'omitnan'),'gaussian',75)),'k')
    %ylim(cx1,[-0.92,0.55])
    xlim(cx4,[min(jumpTime),max(jumpTime)])
    ylabel(cx4,'bump amp')
    xlabel(cx4,'Time (s)')
    ylim(cx4,[0.6,1.8])
    
    %xlim(cx1,[min(jumpTime),max(jumpTime)])
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',5) + SEM_vf(keepIndex)]; 
    SEMlow = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',5) - SEM_vf(keepIndex)];
    plot(cx2,jumpTime,smoothdata(mean(vfSum,'omitnan'),'gaussian',5),'k','LineWidth',1.5)
    patch(cx2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(cx2,mean(smoothdata(mean(vfSum,'omitnan'),'gaussian',5)),'k')
    ylim(cx2,[0,4])
    xlim(cx2,[min(jumpTime),max(jumpTime)])
    ylabel(cx2,'forward velocity (mm/s)')
    xlabel(cx2,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',75) + SEM_vy(keepIndex)]; 
    SEMlow = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vy(keepIndex)];
    plot(cx3,jumpTime,smoothdata(mean(vySum,'omitnan'),'gaussian',75),'k','LineWidth',1.5)
    patch(cx3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(cx3,0,'k')
    ylim(cx3,[-60,60])
    xlim(cx3,[min(jumpTime),max(jumpTime)])
    ylabel(cx3,'rotational velocity (deg/s)')
    xlabel(cx3,'Time (s)')
    
%     keepIndex = ~isnan(SEM_vs);
%     SEMhigh = [smoothdata(mean(vsSum(:,keepIndex),'omitnan'),'gaussian',75) + SEM_vs(keepIndex)]; 
%     SEMlow = [smoothdata(mean(vsSum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vs(keepIndex)];
%     plot(cx4,jumpTime,smoothdata(mean(vsSum,'omitnan'),'gaussian',75),'k','LineWidth',1.5)
%     patch(cx4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     yline(cx4,0,'k')   
%     ylim(cx4,[-1,1])
%     xlim(cx4,[min(jumpTime),max(jumpTime)])
     linkaxes([cx1,cx2,cx3,cx4],'x')

end