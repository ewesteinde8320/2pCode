function PFL3_RLDiff_rho_plots(rootDir)
    foldersStruct = get_folders(rootDir,0,1);

    folders = string({foldersStruct.folder})'; 
    folders = folders(contains(folders,'LAL')); 

    headers = {'folder','maxDiff','rho','estGoal','estAntiGoal','theta'};
    pattern_summary = cell2table(cell(1,6),'VariableNames',headers); 
    cnt = 1; 
    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
    countFail = 1;
    for ff = 1:folderNum
        try
          %% Get folder information
          folder = folders(ff);
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end

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
            clear ftT ZData 
            data_filelist = dir(processedData_dir);
            for files = 1:length(data_filelist)
                if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                    load(fullfile(processedData_dir,data_filelist(files).name));
                end
                load(fullfile(processedData_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']))
            end

            [~, ~, rho_all, theta_all] = CalculateAverageHeading_wholeSegment(ftT,30,threshold,'all', 60);
            theta_all = -theta_all; 
            %%
            total_mov_mm = abs(ftT.velFor{1}) + abs(ftT.velSide{1}) + abs(ftT.velYaw{1}*4.5);
            
            if sum(total_mov_mm > 3)/60 > 60 % moved for at least 60 seconds
                no0vel_idx = find(total_mov_mm > threshold);
                angle = -ftT.cueAngle{1}(no0vel_idx);

                window = 20; 
                edges_angle = [-180 - window/2:window:180 + window/2]; 
                
                activity = ZData.(3){2}(no0vel_idx) - ZData.(3){1}(no0vel_idx);  

                if isnan(activity) | isnan(angle)
                    activity = activity( ~isnan(activity) & ~isnan(angle)); 
                    angle = angle( ~isnan(activity) & ~isnan(angle)); 
                end
                
                [zscore, centers_angle] = binData(activity, angle, edges_angle);
                roiActivity = zscore; 
                %roiActivity = activity; 

                roiActivity(1) = (roiActivity(1) + roiActivity(end))./2;
                roiActivity(end) = []; 
                centers_angle(end) = [];
                
                smoothZ = smoothdata(roiActivity,'gaussian',5);
                figure(33);clf;plot(centers_angle,roiActivity);hold on;yline(0)
                centers_angle(1) = abs(centers_angle(1));
                [bump_params, ~, ~] = fit_sinusoid(roiActivity,[0,2*pi], 1);
                
%                 %% find transition points, where R = L
%                 
                 tDiff = []; 
%                 count = 1; 
%                 for a = 1:length(smoothZ)
%                     if a == 1
%                         prevDiff = smoothZ(end); 
%                         nextDiff = smoothZ(a+1);
%                     elseif a == length(smoothZ)
%                         prevDiff = smoothZ(a-1);
%                         nextDiff = smoothZ(1); 
%                     else
%                         prevDiff = smoothZ(a-1);
%                         nextDiff = smoothZ(a+1);
%                     end   
%                     
%                     if (nextDiff == abs(nextDiff) && prevDiff == abs(prevDiff)) || (nextDiff ~= abs(nextDiff) && prevDiff ~= abs(prevDiff))
%                         tDiff(count) = 0;
%                     else
%                         tDiff(count) = 1; 
%                     end                    
%                     count = count + 1;                   
%                 end
%                         
%                % clean up transition points
%                for t = 1:length(tDiff)
%                    if tDiff(t) == 1 && tDiff(t+1) == 1
%                       [~,minIdx] = min([abs(smoothZ(t)),abs(smoothZ(t+1))]);
%                       if minIdx == 1
%                           tDiff(t+1) = 0;
%                       else
%                           tDiff(t) = 0; 
%                       end
%                    end
%                end
%                
%                % identify which transition points would create fixed points
%                count = 1; 
%                count2 = 1; 
%                fixedPointsIdx = [];
%                unstablePointsIdx = [];
% %                
%                 for a = 1:length(tDiff)
%                    if tDiff(a) == 1
%                         if a == 1
%                             prevDiff = smoothZ(end); 
%                             nextDiff = smoothZ(a+1);
%                         elseif a == length(smoothZ)
%                             prevDiff = smoothZ(a-1);
%                             nextDiff = smoothZ(1); 
%                         else
%                             prevDiff = smoothZ(a-1);
%                             nextDiff = smoothZ(a+1);
%                         end 
%                         
%                         if prevDiff > 0 && nextDiff < 0 
%                             fixedPointsIdx(count) = a;
%                             count = count + 1;
%                         elseif prevDiff < 0 && nextDiff > 0
%                             unstablePointsIdx(count2) = a;
%                             count2 = count2 + 1; 
%                         end
%                    end
%                 end
%                 estGoal = [];
%                 estAntiGoal = []; 
%                 if length(fixedPointsIdx) > 1
%                     % figure out how to handle this
%                     disp('agh')
%                     bGoalArray = ones(size(fixedPointsIdx)) * theta_all;                    
%                     goalDiff = angdiff(bGoalArray,deg2rad(fixedPointsIdx));   
%                     
%                     estGoal = centers_angle(fixedPointsIdx(goalDiff == min(goalDiff))); 
%                     estAntiGoal = wrapTo180(estGoal + 180);
%                     
%                     
%                     estGoal_array = ones(size(centers_angle))*estGoal;
%                     estAntiGoal_array = ones(size(centers_angle))*estAntiGoal;
% 
%                     goalDiff = rad2deg(angdiff(deg2rad(estGoal_array),deg2rad(centers_angle))); % - = counterclockwise to goal, + = clockwise to goal
% 
%                     maxCounterDiff = max(roiActivity(goalDiff < 0)); 
%                     maxClockwiseDiff = abs(min(roiActivity(goalDiff > 0 & goalDiff < 180))); 
% 
%                     maxDiff = maxCounterDiff + maxClockwiseDiff; 
%                     %maxDiff = max(roiActivity) - min(roiActivity);
%                     %maxDiff = max(smoothZ) - min(smoothZ);
% 
%                     pattern_summary.folder(cnt) = {folder};
%                     pattern_summary.maxDiff(cnt) = {maxDiff};
%                     pattern_summary.rho(cnt) = {rho_all}; 
%                     pattern_summary.estGoal(cnt) = {estGoal};
%                     pattern_summary.estAntiGoal(cnt) = {estAntiGoal};
%                     pattern_summary.theta(cnt) = {rad2deg(theta_all)}; 
% 
%                     cnt = cnt + 1; 
%                 elseif length(fixedPointsIdx) == 1
%                     estGoal = centers_angle(fixedPointsIdx); 
%                     estAntiGoal = wrapTo180(estGoal + 180);
%                     
%                     
%                     estGoal_array = ones(size(centers_angle))*estGoal;
%                     estAntiGoal_array = ones(size(centers_angle))*estAntiGoal;
% 
%                     goalDiff = rad2deg(angdiff(deg2rad(estGoal_array),deg2rad(centers_angle))); % - = counterclockwise to goal, + = clockwise to goal
% 
%                     maxCounterDiff = max(roiActivity(goalDiff < 0)); 
%                     maxClockwiseDiff = abs(min(roiActivity(goalDiff > 0 & goalDiff < 180))); 

                    %maxDiff = maxCounterDiff + maxClockwiseDiff; 
                    maxDiff = max(roiActivity) - min(roiActivity);
                    %maxDiff = max(smoothZ) - min(smoothZ);

                    pattern_summary.folder(cnt) = {folder};
                    pattern_summary.maxDiff(cnt) = {maxDiff};
                    pattern_summary.bumpAmp(cnt) = {bump_params.amp};
                    pattern_summary.rho(cnt) = {rho_all}; 
%                     pattern_summary.estGoal(cnt) = {estGoal};
%                     pattern_summary.estAntiGoal(cnt) = {estAntiGoal};
                    pattern_summary.theta(cnt) = {rad2deg(theta_all)}; 

                    cnt = cnt + 1; 
                %end 
            end
        catch
             disp([folder, ' failed'])
        end  
    end
        
        
%         goalNum = [];
%         for f = 1:length(pattern_summary.estGoal)
%             if isempty(pattern_summary.estGoal{f}) 
%                 goalNum(f) = 0;
%             elseif length(pattern_summary.estGoal{f}) == 1
%                 goalNum(f) = 1;
%             else
%                 goalNum(f) = 2; 
%             end
%         end
        
        figure();
        set(gcf,'color','w','renderer','painters')
        hold on
        scatter(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.maxDiff),'k')
%         scatter(cell2mat(pattern_summary.rho(goalNum == 0)),cell2mat(pattern_summary.maxDiff(goalNum == 0)),'k')
%         scatter(cell2mat(pattern_summary.rho(goalNum == 1)),cell2mat(pattern_summary.maxDiff(goalNum == 1)),'b')
%         scatter(cell2mat(pattern_summary.rho(goalNum == 2)),cell2mat(pattern_summary.maxDiff(goalNum == 2)),'r')
        [p,S] = polyfit(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.maxDiff),1); 
        x = linspace(0,max(cell2mat(pattern_summary.rho)),1000);
        fit1 = polyval(p,x);
        plot(x,fit1)
        [R,P] = corrcoef(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.maxDiff));
        annotation('textbox',[0.5,0.8,0.1,0.1],'String',['R = ',num2str(R(1,2)),' p = ',num2str(P(1,2))], 'EdgeColor','none')
        xlabel('rho')
        ylabel('R-L PFL3')
        
        figure();
        set(gcf,'color','w','renderer','painters')
        hold on
        scatter(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.bumpAmp))
        [p,S] = polyfit(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.bumpAmp),1); 
        x = linspace(0,max(cell2mat(pattern_summary.rho)),1000);
        fit1 = polyval(p,x);
        plot(x,fit1)
        [R,P] = corrcoef(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.bumpAmp));
        annotation('textbox',[0.5,0.8,0.1,0.1],'String',['R = ',num2str(R(1,2)),' p = ',num2str(P(1,2))], 'EdgeColor','none')
        xlabel('rho')
        ylabel('bumpAmp')
        
%         figure();
%         scatter(cell2mat(pattern_summary.rho),abs(cell2mat(pattern_summary.angDiff_Z)))
%         
%         hold on
%         scatter(cell2mat(pattern_summary.rho),abs(cell2mat(pattern_summary.angDiff_bump)))

    %     angleShift = wrapTo180(wrapTo180(angle)-wrapTo180(primeAngle));
    % 
    %     window = 90; 
    %     edges_angle = [-180 - window/2:window:180 + window/2]; 
    % 
    %     roiActivity = []; 
    %     for roi = 1:size(Z,2)
    %         activity = Z(:,roi);
    %         if isnan(activity) | isnan(angleShift)
    %             activity = activity( ~isnan(activity) & ~isnan(angleShift)); 
    %             angleShift = angleShift( ~isnan(activity) & ~isnan(angleShift)); 
    %         end
    %         [zscore, centers_angle] = binData(activity, angleShift, edges_angle);
    %         roiActivity(roi,:) = zscore; 
    %     end
    % 
    %     roiActivity(:,1) = (roiActivity(:,1) + roiActivity(:,end))./2;
    %     roiActivity(:,end) = []; 
    % 
    %     activityTable = Z;
    %     ncol = size(roiActivity,2);
    %     color = (cbrewer2('Greens', ncol));
    % 
    %     test = circshift(color,8,1);
    % 
    %     brainSpace = linspace(180,-180,size(Z,2));
    %     legendNames = []; 
    %     legendNames(:,1)= centers_angle(1:end-1);
    %     legendNames = num2str(legendNames);
    % 
    %     figure()
    %     set(gcf,'color','w','renderer','painters')
    %     hold on
    %     for h = 1:size(roiActivity,2)
    %         plot(brainSpace,roiActivity(:,h),'color',test(h ,:),'LineWidth',1.5)
    %     end
    %     legend({legendNames})
    %     xlabel('brain space')
    %     ylabel('z-scored df/f')
    %     legend boxoff
    %     title(num2string(rho_all))


end





