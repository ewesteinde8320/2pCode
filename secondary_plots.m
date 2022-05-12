function secondary_plots(rootDir, savePlots)

%% Default settings

    arguments
        rootDir char
        savePlots logical = 1 
    end

    
%% get folder list

    folders = get_folders(rootDir);
           
    
    %% Process each folder
    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
    for ff = 1:folderNum
      
      %% Get folder information
      folder = folders(ff).folder;
      
       
%% Load in fictrac & ROI data
        try
            if strcmp(folder(end),'.')
                folder = folder(1:end-2); 
            end
            
            %imaging_fictrac_pipeline(folder, 0, 0, 0, 1);
            
            % Make images directory
            imagesDir = fullfile(folder,'images');
            if ~exist(imagesDir, 'dir')
                mkdir(imagesDir)
            end

            lineplotDir = fullfile(imagesDir,'lineplots');
            if ~exist(lineplotDir, 'dir')
                mkdir(lineplotDir)
            end
            
            heatmapDir = fullfile(imagesDir,'heatmaps');
            if ~exist(heatmapDir, 'dir')
                mkdir(heatmapDir)
            end

             processedData_dir = fullfile(folder,'processed_data');

            if ~exist(processedData_dir,'dir')
                disp('No processed data found, running process2p_fictrac_data now')
                process2p_fictrac_data(folder)
            end

            % Get data files
            expID = get_expID(folder);
            expList = {expID};

            % Load metadata 
            [expMd, trialMd] = load_metadata(expList, folder);

            % Load imaging data
            roiData = load_roi_data(expList, folder);

            % Load FicTrac data
            ftData = load_ft_data(expList, folder);

            % Load panels metadata
            panelsMetadata = load_panels_metadata(expList, folder);


            numTrials = max(size(unique(roiData.trialNum),1),length(trialMd.trialNum)); 
        %%
        for nTrial = 1:numTrials


         %% extract fictrac Data
         % Get non downsmapled fictrac values for heading dist plots

        if any(ftData.trialNum==nTrial)
            ftT = ftData(ftData.trialNum ==nTrial, :);
            ftT_time = ftT.frameTimes{1};

        %% Plot heading 
            ftT_cueAngle = ftT.cueAngle{1};
            ftT_cueAngle_smooth = ftT_cueAngle;

            %% Plot fwSpeed 
            ftT_fwSpeed = ftT.fwSpeed{1};
            %ftT_fwSpeed_smooth = smoothdata(ftT_fwSpeed, 1, 'loess', 10);

            %% Plot sideSpeed 
            ftT_sideSpeed = ftT.sideSpeed{1};
            %ftT_sideSpeed_smooth = smoothdata(ftT_sideSpeed, 1, 'loess', 10);

            %% Plot yawSpeed 
            ftT_yawSpeed = ftT.yawSpeed{1};
            %ftT_yawSpeed_smooth = smoothdata(ftT_yawSpeed, 1, 'loess', 10);
        end

        %% get processed data

            data_filelist = dir(processedData_dir);
            for files = 1:length(data_filelist)
                if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                    load(fullfile(processedData_dir,data_filelist(files).name));
                end
            end

            roiNames = roiData.roiName;

        %% plot activty vs behaviour
        
        Activity_Heatmap(ftT_down, roiData, Z, Zf, nTrial,savePlots, heatmapDir, expID)

        activityVSbehaviour_lineplots(ftT_down, Z, Zf, roiData, nTrial, expMd, savePlots,lineplotDir, expID)
        threshold = 3;
        activityVSbehaviour_no0vel_secplots(ftT_down, Z, Zf, roiData, nTrial, threshold, expMd,savePlots,lineplotDir, expID)

        %% heading distribution plots
        behaviourData.vel_for = ftT_fwSpeed; 
        behaviourData.vel_side = ftT_sideSpeed;
        behaviourData.vel_yaw = wrapTo180((ftT_yawSpeed/ (2*pi) ) * 360);
        behaviourData.angle = ftT_cueAngle_smooth;
        minVel = 1.5; 
        window = 60; 
        sampRate = 60; %fictrac data samp rate 

        [rho, theta, ~, f, g] = plotHeadingDist(window, minVel, behaviourData, sampRate, 1);


        if savePlots == 1
            saveas(f, fullfile(imagesDir,[expID,'_',num2str(nTrial),'_headingVectors.fig']));
            saveas(g, fullfile(imagesDir,[expID,'_',num2str(nTrial),'_headingDist.fig']));
        %clear
        end
        measurement = 1;
        interactive = 0; 

        %% replot activity-behaviour relationship at diff headings, must be run after heading dis plots
        activityVSbehaviour_headings(rho, theta, roiData, Z, Zf, ftT_down, trialMd, expMd, nTrial, measurement, savePlots, interactive,lineplotDir, expID)  
        activityVSbehaviour_headings_no0vel(rho, theta, roiData, threshold, Z, Zf,ftT_down, trialMd, expMd, nTrial, measurement, savePlots, interactive,lineplotDir,expID)  
        
        if regexp(folder, 'LAL')
            activityVSbehaviour_headings_sepROIs(rho, theta, roiData, Z, ftT_down, nTrial, measurement, savePlots,interactive,lineplotDir, expID)
        end
        
        % %% calculate bump position & offset
        % in progress, not working yet
        % midline distances = vector of midline distances
        % dff_data = array of dff values where row = timepoint col = midline
        % distance

        % 
        % if ~contains(expMd.expName{1},'LAL')
        %     if contains(expMd.expName{1},'FB')
        %         roi_root = 'FB_row_line_';
        %     elseif contains(expMd.expName{1},'PB')
        %         roi_root = 'PB_row_line_';
        %     else
        %         error('unknown brain region')
        %     end
        %         [~,midline_distances_idx] = regexp(org_roiData.roiName,roi_root);
        %         midline_distances_idx = cell2mat(midline_distances_idx);
        %         roiNames = org_roiData.roiName{:};
        %         for name = 1:size(org_roiData.roiName,1)
        %             midline_distances(name) = str2double(org_roiData.roiName{name}(midline_distances_idx + 1:end));
        %             dff_data(:,name) = flData_all{name}; 
        %         end
        %     
        % 
        %     [bump_mag,bump_width,adj_rs,bump_pos] = fitVonMises(midline_distances,dff_data);
        %end
        end
        pause(2)
        close all
        catch
            disp([folder,' failed to process'])
        end
    end
end