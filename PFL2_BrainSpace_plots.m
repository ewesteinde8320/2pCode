function PFL2_BrainSpace_plots(rootDir)
    foldersStruct = get_folders(rootDir,1,0);

    folders = string({foldersStruct.folder})'; 
    folders = [folders(contains(folders,'PB'));folders(contains(folders,'FB'))]; 

    headers = {'folder','angDiff','ampDiff_bump','ampVar_bump','ampDiff_z','ampVar_z','rho','estGoal','estAntiGoal','theta'};
    pattern_summary = cell2table(cell(1,10),'VariableNames',headers); 
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

            data_filelist = dir(processedData_dir);
            for files = 1:length(data_filelist)
                if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                    load(fullfile(processedData_dir,data_filelist(files).name));
                end
                load(fullfile(processedData_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']))
            end

            [~, ~, rho_all, theta_all] = CalculateAverageHeading_wholeSegment(ftT,30,threshold,'all', 60);
            %%
            total_mov_mm = abs(ftT.velFor{1}) + abs(ftT.velSide{1}) + abs(ftT.velYaw{1}*4.5);
            
            if sum(total_mov_mm > threshold)/60 > 60 % moved for at least 60 seconds
                no0vel_idx = find(total_mov_mm > threshold);
                angle = -ftT.cueAngle{1}(no0vel_idx);

                window = 2; 
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
                centers_angle(:,1) = abs(centers_angle(:,1));
                centers_angle(:,end) = [];

                [bump_params, ~, ~] = fit_sinusoid(roiActivity,[0,2*pi], 0);

%                 estGoal = centers_angle(bump_params.amp == min(bump_params.amp));
%                 estAntiGoal = centers_angle(bump_params.amp == max(bump_params.amp)); 
                
                for b = 1:length(centers_angle)/2
                    distances(b) = diff([bump_params.amp(b), bump_params.amp(b+length(centers_angle)/2)]); 
                end
                
                maxDiffidx = find(abs(distances) == max(abs(distances)));
                
                if distances(maxDiffidx) > 0 
                    estGoal = centers_angle(maxDiffidx);
                    estAntiGoal = centers_angle(maxDiffidx + length(centers_angle)/2);
                else
                    estAntiGoal = centers_angle(maxDiffidx);
                    estGoal = centers_angle(maxDiffidx + length(centers_angle)/2);
                end
                    
                    
                    
                
                
                %%
                angle = -ftT.cueAngle{1}(no0vel_idx);
                angle = wrapTo180(wrapTo180(angle)-wrapTo180(estGoal));
                window = 90; 
                edges_angle = [-180 - window/2:window:180 + window/2];
                
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
                centers_angle(:,1) = abs(centers_angle(:,1));
                centers_angle(:,end) = [];

                [bump_params, ~, ~] = fit_sinusoid(roiActivity,[0,2*pi], 0);

                 
                %angDiff_bump = rad2deg(angdiff(deg2rad(centers_angle(bump_params.amp == max(bump_params.amp))),deg2rad(centers_angle(bump_params.amp == min(bump_params.amp)))));
                ampDiff_bump = max(bump_params.amp(1)) - min(bump_params.amp(3));
                Z_ranges = max(roiActivity,[],1)-min(roiActivity,[],1);
                ampDiff_Z = Z_ranges(1) - Z_ranges(3); 
                %angDiff_Z = rad2deg(angdiff(deg2rad(centers_angle(Z_ranges == max(Z_ranges))),deg2rad(centers_angle(Z_ranges == min(Z_ranges)))));
                ampVar = var(bump_params.amp); 


                pattern_summary.folder(cnt) = {folder};
                %pattern_summary.angDiff_bump(cnt) = {angDiff_bump};
                pattern_summary.ampDiff_bump(cnt) = {ampDiff_bump};
                %pattern_summary.angDiff_Z(cnt) = {angDiff_Z};
                pattern_summary.ampDiff_Z(cnt) = {ampDiff_Z};
                pattern_summary.ampVar(cnt) = {ampVar};
                pattern_summary.rho(cnt) = {rho_all}; 
                pattern_summary.estGoal(cnt) = {estGoal};
                pattern_summary.estAntiGoal(cnt) = {estAntiGoal};
                pattern_summary.theta(cnt) = {-rad2deg(theta_all)}; % rho is in cue pos coord not HD, correct for this

                cnt = cnt + 1; 

                %% Plot activity pattern across brain space, 0 = HD with highest amp sinusoid across PFL2 neurons
                disp(folder)
                %plot_brainSpace_activity_byHeading(folder, estGoal,rho_all)

                %pause(0.002)
            end
        catch
             disp([folder, ' failed'])
        end  
    end
         for f = 1:size(pattern_summary,1)
        end
        
        figure();
        set(gcf,'color','w','renderer','painters')
        hold on
        scatter(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_bump))
        [p,S] = polyfit(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_bump),1); 
        x = linspace(0,max(cell2mat(pattern_summary.rho)),1000);
        fit1 = polyval(p,x);
        plot(x,fit1)
        [R,P] = corrcoef(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_bump));
        annotation('textbox',[0.5,0.8,0.1,0.1],'String',['R = ',num2str(R(1,2)),' p = ',num2str(P(1,2))], 'EdgeColor','none')
        xlabel('rho')
        ylabel('max - min bump amp')
        
        figure();
        set(gcf,'color','w','renderer','painters')
        hold on
        scatter(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_Z))
        [p,S] = polyfit(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_Z),1); 
        x = linspace(0,max(cell2mat(pattern_summary.rho)),1000);
        fit1 = polyval(p,x);
        plot(x,fit1)
        [R,P] = corrcoef(cell2mat(pattern_summary.rho),cell2mat(pattern_summary.ampDiff_Z));
        annotation('textbox',[0.5,0.8,0.1,0.1],'String',['R = ',num2str(R(1,2)),' p = ',num2str(P(1,2))], 'EdgeColor','none')
        xlabel('rho')
        ylabel('max - min Z')
        
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





