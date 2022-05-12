function process2p_fictrac_data(rootDir)

    filelist = dir(fullfile(rootDir, '**/*.*'));  % get list of files and folders in any subfolder
    filelist = filelist([filelist.isdir]);  % only folders from list
    filelistNum = length(filelist);
    correct = [];
    for f = 1:filelistNum
      baseFolderName = filelist(f).name;
      if strcmp(baseFolderName,'..') || ~strcmp(baseFolderName,'.')
          present = false;
      else
          parentFolderName = filelist(f).folder;
          fullFolderName = fullfile(parentFolderName, baseFolderName);
          filelist(f).folder = fullFolderName;
          contents = fullfile(fullFolderName,'config.txt'); % fictrac data and ...
          present = exist(contents,'file')>0;
          listingTiff = dir(fullfile(fullFolderName,'*_daqData_*.mat'));
          present = numel(listingTiff) > 0 && present; % Image data
      end
      correct(end+1) = present;
    end 
    folders = filelist(logical(correct));
   
    %% Process each folder
    folderNum = length(folders);
    fprintf(1, '##### Found %d potential experiment folders to process...#####\n', folderNum);
    for ff = 1:folderNum
        close all 
        %clearvars -except ff folderNum folders savePlots
      
      %% Get folder information
      parentDir = folders(ff).folder;
      
      if strcmp(parentDir(end),'.') 
          parentDir = parentDir(1:end-2);
      end

%% Load in fictrac & ROI data
 
    % Analysis settings
    p = [];
    p.smWin = 5;
    p.flType = 'expDff';
    ztab = [];
    
    % Get data files
    expID = get_expID(parentDir);
    expList = {expID};
    
    % Load metadata 
    [expMd, trialMd, ~, ~] = load_metadata(expList, parentDir);

    % Load imaging data
    roiData_all = load_roi_data(expList, parentDir);

    % Load FicTrac data
    ftData = load_ft_data(expList, parentDir);
    
    numTrials = max(size(unique(roiData_all.trialNum),1),length(trialMd.trialNum)); 
    

%%
        for nTrial = 1:numTrials
            if any(roiData_all.trialNum==nTrial)

            roiData = roiData_all(roiData_all.trialNum==nTrial,:);
            
            %% Volume rate
            if size(trialMd,1) >= nTrial
                tMd = trialMd(trialMd.trialNum == nTrial, :);
            else 
                tMd = trialMd;
            end 
            
            %% resample ROI data to integer volume rate
            try
                roi_time = [0:1/trialMd.volumeRate:trialMd.trialDuration]';
            catch
                roi_time = [0:1/expMd.volumeRate:trialMd.trialDuration]';
                tMd.volumeRate = expMd.volumeRate;
            end
                
            roi_time = seconds(roi_time);
            
                if length(roi_time) == length(roiData.rawFl{1}) -1 
                    roi_time(end+1) = seconds(600); 
                end
            
            for r = 1:size(roiData,1)
                rawFl = roiData.rawFl{r};
                roi = table(rawFl); 
                roi.seconds = roi_time; 
                roi_timetable_temp = table2timetable(roi,'RowTimes','seconds'); 
                roi_timetable = retime(roi_timetable_temp,'regular','nearest','SampleRate',round(tMd.volumeRate)); 
                roi = timetable2table(roi_timetable(1:end-1,:)); 
                roiData.rawFl{r} = roi.rawFl;
                roiData.time{r}= roi.seconds;
                roiData.sampRate(r) = round(tMd.volumeRate);
            end

            volRate =roiData.sampRate(1);
            ztab = [];
            zftab = []; 

            roiNames = roiData.roiName;
            %% For each ROI plot with the correct ROI number
            count = 0; 
            for nRoi = 1:numel(roiNames)

                %% Get ROI name
                roiName = char(roiNames(nRoi));
                %% If the ROi exists with data 
                currRoiData = roiData(roiData.trialNum == nTrial & strcmp(roiData.roiName, roiName), :);
                if ismember(roiName, roiData.roiName) && ~isempty(currRoiData) && ~contains(currRoiData.roiName,'mid') && ~isempty(currRoiData.rawFl{1}) %temp addition 
                    %% Get dF/F % not downsampled need to determine if this is matched accurately 
                    count = count + 1;
                    flData = deltaFoverF(roiData, nTrial, roiName, p);
                    mflData = Mad(roiData, nTrial, roiName, p);
                    for ii = 1:size(mflData,1)
                        zcell = {count,ii/volRate,mflData(ii)}; %{roiName,ii,mflDataDown(ii)};
                        ztab = [ztab;zcell];
                        zfcell = {count,ii/volRate,flData(ii)};
                        zftab = [zftab;zfcell];
                    end
        %                     else %% But if not ...
        %                         flData = zeros(size(roiTime));
                end
                flData_all{count} = flData;
                flData_names{count} = roiName;
            end
        end

        Z = cell2table(ztab, 'VariableNames',{'roiName' 'second' 'data'});
        Zf = cell2table(zftab, 'VariableNames',{'roiName' 'second' 'data'});

        % match fictrac data to imaging data
        %%
        % strategy 1: bin behaviour data leading up to every volume timepoint &
        % average
%         vf_bin = []; 
%         count = 0; 
%         prev_vol_time = 0; 
%         for vol_time = min(Z.second):1/volRate:max(Z.second)
%             count = count + 1; 
%             f_idx = find(ftT_time > prev_vol_time & ftT_time <= vol_time & ftT_time ~= 0);
%             vf_bin(count) = median(ftT_fwSpeed_smooth(f_idx)); 
%             prev_vol_time = vol_time; 
%         end
          

        %%
        % strategy 2: find closest fictrac timepoint to each vol timepoint
        
        if any(ftData.trialNum==nTrial)
            ftT = ftData(ftData.trialNum ==nTrial, :);
            ftT_time = ftT.frameTimes{1};
        end

%         A = repmat(ftT_time,[1 length(Z.second(Z.roiName == 1))]);
%         [minValue,closestIndex] = min(abs(A-(Z.second(Z.roiName == 1))'));

        %ftT_down = cell2table(cell(0,length(var_names)), 'VariableNames', var_names);
        
%         for col = 1:width(ftT)
%             if iscell(ftT.(col)) && length(ftT.(col){1}) > length(closestIndex)
%                  if cell2mat(regexp(var_names(col),'Speed'))
%                      smooth_table = array2table({smoothdata(ftT.(col){1}, 1, 'loess', 20)});
%                      temp_table = array2table({smooth_table.(1){1}(closestIndex)}); 
%                  else
%                     temp_table = array2table({ftT.(col){1}(closestIndex)}); 
%                  end 
%             else
%                 temp_table = ftT(:,col); 
%             end          
%             ftT_down(:,col) = temp_table;
%         end

        var_names = ftT.Properties.VariableNames;         
        for col = 1:width(ftT)
            if strcmp(var_names{col},'expID')
                temp_table = ftT(:,col); 
            elseif iscell(ftT.(col)) && ~isempty(ftT.(col){1})
                 if cell2mat(regexp(var_names(col),'Speed'))
                     %smooth_table = array2table({smoothdata(ftT.(col){1}, 1, 'loess', 500)});
                     if isnan(ftT.(col){1}(1)) || isinf(ftT.(col){1}(1))
                        ftT.(col){1}(1) = 0;
                     end
                     temp_table = array2table({smoothdata(resample_with_padding(ftT.(col){1},volRate,60)', 1, 'loess', 5)}); 
                 else
                    temp_table = array2table({resample_with_padding(ftT.(col){1},volRate,60)'}); 
                 end 
            else
                temp_table = ftT(:,col); 
            end          
            ftT_down(:,col) = temp_table;
        end

        allVars = 1:width(ftT_down); 
        ftT_down = renamevars(ftT_down,allVars,var_names);
 
%         
%         %% smoothing comparison plots
%         
%         figure()
%         plot(ftData.trialTime{1},ftData.(col){1})
%         hold on
%         plot(ftT_down.trialTime{1},ftT_down.(col){1})
%         hold on
%         plot(ftT_down.trialTime{1},temp_table.(1){1})
%         legend('original','downSampled')
        

        processed_data_dir = fullfile(parentDir,'processed_data');
        if ~exist(processed_data_dir, 'dir')
            mkdir(processed_data_dir)
        end


        save(fullfile(processed_data_dir,['fictracData_down_Trial_00',num2str(nTrial),'.mat']), 'ftT_down','-v7.3'); 
        save(fullfile(processed_data_dir,['df_f_Trial_00',num2str(nTrial),'.mat']), 'Zf','-v7.3');
        save(fullfile(processed_data_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']), 'Z','-v7.3');
        save(fullfile(processed_data_dir,['roiData_Trial_00',num2str(nTrial),'.mat']), 'roiData','-v7.3');

        ftCSV = flyg2csv(ftT_down,'intX');
        fileCSV = fullfile(processed_data_dir,['ftData_down_Trial_00',num2str(nTrial),'.csv']);
        writetable(ftCSV,fileCSV,'Delimiter',',','QuoteStrings',true);
        disp(['Wrote: ' fileCSV]);
    
        ftCSV = flyg2csv(roiData,'rawFl');
        fileCSV = fullfile(processed_data_dir,['roiData_Trial_00',num2str(nTrial),'.csv']);
        writetable(ftCSV,fileCSV,'Delimiter',',','QuoteStrings',true);
        disp(['Wrote: ' fileCSV]);

        end
    end
    
    
    
    % testing plot
%     figure(2);clf;
%     plot(ftT.frameTimes{1}, ftT.fwSpeed{1});
%     hold on
%     plot(ftT_down.frameTimes{1},ftT_down.fwSpeed{1}); 
%     legend
    
    disp('-------Finished processing Data----------')
end    
    
    

    
    
    
    
    