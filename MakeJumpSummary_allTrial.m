function MakeJumpSummary_allTrial(rootDir,PFL2,PFL3)       

     folders = get_folders(rootDir,PFL2,PFL3);
     
     headers = {'folder','trial','jumpSize','jumpIdx','velFor','velSide','velYaw','cueAngle','corrected','timeMov'};
     jump_summary = cell2table(cell(1,10),'VariableNames',headers); 
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
        
            data_filelist = dir(processedData_dir);
            for files = 1:length(data_filelist)
                if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                    load(fullfile(processedData_dir,data_filelist(files).name));
                end
            end

    %% Remove idx where the fly isn't moving
            %ftT = ftT_minSmooth; 
            window = 10; 
            [~, jump_array, ~] = detect_jumps(ftT, window, 5,0);

                count90 = 1;
                count180 = 1;
                jump90_time = []; 
                jump180_time = []; 
                
                for jump = 1:size(jump_array,1)
                    jumpSize = jump_array(jump,4);
                    jumpIdx = [jump_array(jump,1):jump_array(jump,3)]'; 
                    jump_vf = ftT.velFor{1}(jumpIdx);
                    jump_vy = rad2deg(ftT.velYaw{1}(jumpIdx));
                    jump_vs = ftT.velSide{1}(jumpIdx);
                    jump_angle = ftT.cueAngle{1}(jumpIdx);
                    preJumpIdx = [jump_array(jump,1):jump_array(jump,2)-1];
                    timeMov = sum((abs(ftT.velFor{1}(preJumpIdx))+ abs(ftT.velSide{1}(preJumpIdx)) + abs(ftT.velYaw{1}(preJumpIdx))*4.5) > 1.5)/60; % seconds

                    [rho, preJumpAngle] = CalculateAverageHeading(ftT,1, preJumpIdx);
                    
                    %preJumpAngle = ftT.cueAngle{1}(preJumpIdx(end));
                    
                    postJumpIdx = jump_array(jump,2)+1;
                    sizePostJump = size(ftT.cueAngle{1}(postJumpIdx:jump_array(jump,3)));
                    preJump = ones(sizePostJump) * preJumpAngle;
                    cueDiff = angdiff(preJump,deg2rad(ftT.cueAngle{1}(postJumpIdx:jump_array(jump,3)))); 
                    
                    if jumpSize == 180
                        correct = sum(abs(cueDiff) < deg2rad(37.5)); 
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
                    jump_summary_timeMov(count) = {timeMov}; 

                  
                    if  rho > 0.88 && timeMov >= 1  && correct > 0
                        jump_summary.corrected(count) = {1}; 
                    elseif correct > 0 && rho < 0.88
                        jump_summary.corrected(count) = {2}; 
                    elseif correct == 0 && timeMov >= 1
                        jump_summary.corrected(count) = {0}; 
                    else
                        jump_summary.corrected(count) = {3}; 
                    end
                    
                    count = count + 1; 
                    
                end
        end
%         for jump = 1:size(jump_summary,1)
%             figure();
%             plot(jump_summary.cueAngle{jump})
%         end

        catch
            disp([folder, ' failed'])
        end
    end  
    Sum180 = jump_summary(cell2mat(jump_summary.jumpSize) == 180 & cell2mat(jump_summary.corrected) == 1 ,:); 
    Sum90 = jump_summary(abs(cell2mat(jump_summary.jumpSize)) == 90 & cell2mat(jump_summary.corrected) == 1 ,:);
    Summ90 = jump_summary(cell2mat(jump_summary.jumpSize) == -90 & cell2mat(jump_summary.corrected) == 1 ,:);
 %%   
    jumpTime = [-10+(1/60):1/60:10+(1/60)]; 
    vfSum = zeros(size(Sum180,1),size(jump_vf,1));
    vySum = zeros(size(Sum180,1),size(jump_vf,1));
    vsSum = zeros(size(Sum180,1),size(jump_vf,1));
    figure();
   
    sgtitle('180 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
%     ax2 = subplot(2,1,1);
%     hold on
%     ax3 = subplot(2,1,2); 
%     hold on
%     ax4 = subplot(3,1,3); 
%     hold on
    for jump = 1:size(Sum180,1)  
        %plot(ax1, jumpTime,Sum180.cueAngle{jump})
        %plot(ax2,jumpTime,Sum180.velFor{jump})
        vfSum(jump,:) = Sum180.velFor{jump}'; 
        %plot(ax3,jumpTime,abs(Sum180.velYaw{jump}))
        vySum(jump,:) = abs(Sum180.velYaw{jump})' - abs(Sum180.velYaw{jump}(1));
        %plot(ax4,jumpTime,abs(Sum180.velSide{jump}))
        vsSum(jump,:) = abs(Sum180.velSide{jump})';
    end
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vs = std(vsSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    
%     %xlim(ax1,[min(jumpTime),max(jumpTime)])
%     keepIndex = ~isnan(SEM_vf);
%     SEMhigh = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',75) + SEM_vf(keepIndex)]; 
%     SEMlow = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vf(keepIndex)];
%     plot(ax2,jumpTime,smoothdata(mean(vfSum,'omitnan'),'gaussian',75),'k','LineWidth',1.5)
%     patch(ax2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     %yline(ax2,mean(smoothdata(mean(vfSum,'omitnan'),'gaussian',75)),'k')
%     ylim(ax2,[0,3])
%     xlim(ax2,[min(jumpTime),max(jumpTime)])
%     ylabel(ax2,'forward velocity (mm/s)')
%     xlabel(ax2,'Time (s)')
    smoothVy = smoothdata(mean(vySum,'omitnan'),'gaussian',75); 
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',75) + SEM_vy(keepIndex)-smoothVy(1)]; 
    SEMlow = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vy(keepIndex)-smoothVy(1)];
    plot(jumpTime,smoothdata(mean(vySum,'omitnan'),'gaussian',75)-smoothVy(1),'k','LineWidth',1.5)
    patch([jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(ax3,mean(smoothdata(mean(vySum,'omitnan'),'gaussian',75)),'k')
    %ylim([0,80])
    xlim([min(jumpTime),max(jumpTime)])
    ylabel('rotational speed (deg/s)')
    xlabel('Time (s)')
    
%     keepIndex = ~isnan(SEM_vs);
%     SEMhigh = [smoothdata(mean(vsSum(:,keepIndex),'omitnan'),'gaussian',75) + SEM_vs(keepIndex)]; 
%     SEMlow = [smoothdata(mean(vsSum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vs(keepIndex)];
%     plot(ax4,jumpTime,smoothdata(mean(vsSum,'omitnan'),'gaussian',75),'k','LineWidth',1.5)
%     patch(ax4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     yline(ax4,mean(smoothdata(mean(vsSum,'omitnan'),'gaussian',75)),'k')
%     ylim(ax4,[0,1.5])
%     xlim(ax4,[min(jumpTime),max(jumpTime)])
%     linkaxes([ax2,ax3,ax4],'x')
        
    vfSum = zeros(size(Sum90,1),size(jump_vf,1));
    vySum = zeros(size(Sum90,1),size(jump_vf,1));
    vsSum = zeros(size(Sum90,1),size(jump_vf,1));
%     figure();
%     sgtitle('90 jumps')
%     set(gcf,'color','w')
%     set(gcf,'renderer','painters')
%     bx2 = subplot(2,1,1);
%     hold on
%     bx3 = subplot(2,1,2);
%     hold on
%     bx4 = subplot(3,1,3); 
%     hold on
    for jump = 1:size(Sum90,1)  
        %plot(bx1, jumpTime,Sum90.cueAngle{jump})
        %plot(bx2,jumpTime,Sum90.velFor{jump})
        vfSum(jump,:) = Sum90.velFor{jump}'; 
        %plot(bx3,jumpTime,Sum90.velYaw{jump})
        vySum(jump,:) = abs(Sum90.velYaw{jump}');
        %plot(bx4,jumpTime,Sum90.velSide{jump})
        vsSum(jump,:) = Sum90.velSide{jump}';
    end
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vs = std(vsSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
%     keepIndex = ~isnan(SEM_vf);
%     SEMhigh = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',75) + SEM_vf(keepIndex)]; 
%     SEMlow = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vf(keepIndex)];
%     %xlim(bx1,[min(jumpTime),max(jumpTime)])
%     plot(bx2,jumpTime,smoothdata(mean(vfSum,'omitnan'),'gaussian',75),'k','LineWidth',1.5)
%     patch(bx2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     yline(bx2,mean(smoothdata(mean(vfSum,'omitnan'),'gaussian',75)),'k')
%     ylim(bx2,[0,3])
%     xlim(bx2,[min(jumpTime),max(jumpTime)])
%     ylabel(bx2,'forward velocity (mm/s)')
%     xlabel(bx2,'Time (s)')
    
%     baselineVel = mean(vySum(:,1:500),2,'omitnan');
%     baselineVel = vySum(:,1);
    baselineVel = mean(vySum,1,'omitnan');
    baselineVel = baselineVel(1);
    
    hold on
    keepIndex = ~isnan(SEM_vy);
    smoothVy = smoothdata(mean(vySum,'omitnan'),'gaussian',75); 
    SEMhigh = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',75)  + SEM_vy(keepIndex)]-smoothVy(1); 
    SEMlow = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vy(keepIndex)]-smoothVy(1);
    plot(jumpTime,smoothdata(mean(vySum,'omitnan'),'gaussian',75) - smoothVy(1),'b','LineWidth',1.5)
    patch([jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0,0,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    %yline(0,'k')
    xlim([min(jumpTime),max(jumpTime)])
    ylabel('rotational velocity (deg/s)')
    xlabel('Time (s)')
    
%     keepIndex = ~isnan(SEM_vs);
%     SEMhigh = [smoothdata(mean(vsSum(:,keepIndex),'omitnan'),'gaussian',75) + SEM_vs(keepIndex)]; 
%     SEMlow = [smoothdata(mean(vsSum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vs(keepIndex)];
%     plot(bx4,jumpTime,smoothdata(mean(vsSum,'omitnan'),'gaussian',75),'k','LineWidth',1.5)
%     patch(bx4,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
%     yline(bx4,0,'k')
%     ylim(bx4,[-1,1])
%     xlim(bx4,[min(jumpTime),max(jumpTime)])
%     linkaxes([bx2,bx3,bx4],'x')
 box off


    
    


  %%  
    vfSum = zeros(size(Summ90,1),size(jump_vf,1));
    vySum = zeros(size(Summ90,1),size(jump_vf,1));
    vsSum = zeros(size(Summ90,1),size(jump_vf,1)); 
    figure();
    box off
    sgtitle('-90 jumps')
    set(gcf,'color','w')
    set(gcf,'renderer','painters')
    %cx1 = subplot(4,1,1);
    %hold on
    cx2 = subplot(2,1,1);
    hold on
    cx3 = subplot(2,1,2); 
    hold on
%     cx4 = subplot(3,1,3); 
%     hold on
    for jump = 1:size(Summ90,1)  
        %plot(cx1, jumpTime,Summ90.cueAngle{jump})
        %plot(cx2,jumpTime,Summ90.velFor{jump})
        vfSum(jump,:) = Summ90.velFor{jump}'; 
        %plot(cx3,jumpTime,Summ90.velYaw{jump})
        vySum(jump,:) = Summ90.velYaw{jump}';
        %plot(cx4,jumpTime,Summ90.velSide{jump})
        vsSum(jump,:) = Summ90.velSide{jump}';
    end
    SEM_vf = std(vfSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vy = std(vySum,[],1,'omitnan') / sqrt(size(vfSum,1));
    SEM_vs = std(vsSum,[],1,'omitnan') / sqrt(size(vfSum,1));
    
    %xlim(cx1,[min(jumpTime),max(jumpTime)])
    keepIndex = ~isnan(SEM_vf);
    SEMhigh = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',75) + SEM_vf(keepIndex)]; 
    SEMlow = [smoothdata(mean(vfSum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vf(keepIndex)];
    plot(cx2,jumpTime,smoothdata(mean(vfSum,'omitnan'),'gaussian',75),'k','LineWidth',1.5)
    patch(cx2,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    yline(cx2,mean(smoothdata(mean(vfSum,'omitnan'),'gaussian',75)),'k')
    ylim(cx2,[0,3])
    xlim(cx2,[min(jumpTime),max(jumpTime)])
    ylabel(cx2,'forward velocity (mm/s)')
    xlabel(cx2,'Time (s)')
    
    keepIndex = ~isnan(SEM_vy);
    SEMhigh = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',75) + SEM_vy(keepIndex)]; 
    SEMlow = [smoothdata(mean(vySum(:,keepIndex),'omitnan'),'gaussian',75) - SEM_vy(keepIndex)];
    plot(cx3,jumpTime,smoothdata(mean(vySum,'omitnan'),'gaussian',75),'k','LineWidth',1.5)
    patch(cx3,[jumpTime(keepIndex) fliplr(jumpTime(keepIndex))],[SEMhigh fliplr(SEMlow)],[0.75,0.75,0.75],'FaceAlpha',0.25,'EdgeColor','none')
    yline(cx3,0,'k')
    ylim(cx3,[-35,35])
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
%     linkaxes([cx2,cx3,cx4],'x')

end