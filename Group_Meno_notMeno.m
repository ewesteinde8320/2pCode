function [MenoData, nMenoData] = Group_Meno_notMeno(rootDir)
varNames = ["Folder","numTrial","Goal","rho","timeMov","Indices"];
varTypes = ["string","double","double","double","double","cell"];
MenoData = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);
nMenoData = table('Size',[1 6],'VariableNames', varNames,'VariableTypes',varTypes);

PFL3 = 0; 
PFL2 = 1; 

folders = get_folders(rootDir, PFL2, PFL3);

window = 30; 
minVel = 1.5; 
highThres = 0.88; 
lowThres = 0.88; 

mCount = 1; 
nCount = 1; 

for f = 1:size(folders,1)
    folder = folders(f).folder; 
    if strcmp(folder(end),'.')
        folder = folder(1:end-2); 
    end
    
    expID = get_expID(folder);
    expList = {expID};
    
    [~,ftData, ~] = load_ft_data(expList, folder, 1, 0);

    % Load metadata 
    [expMd, trialMd] = load_metadata(expList, folder);

    % Load imaging data
    roiData = load_roi_data(expList, folder);
    
    processedData_dir = fullfile(folder,'processed_data');
    try
        numTrials = max(size(unique(roiData.trialNum),1),length(trialMd.trialNum));
    catch
        numTrials = 1; 
    end
    
    for nTrial = 1:numTrials
        
        data_filelist = dir(processedData_dir);
        for files = 1:length(data_filelist)
            if regexp(data_filelist(files).name,'.mat') & regexp(data_filelist(files).name,['00',num2str(nTrial)])
                load(fullfile(processedData_dir,data_filelist(files).name));
            end
        end

%%  
    Meno_chunks = [];
    not_Meno_chunks = []; 
    [~, ~, Meno_chunks, not_Meno_chunks] = SegmentMenovsNotMeno(ftT, window, minVel,highThres,lowThres,folder,nTrial,0);
    
    if ~isempty(Meno_chunks) || ~isempty(not_Meno_chunks)
        Meno_chunks = Meno_chunks';
        not_Meno_chunks = not_Meno_chunks';

        count = 1; 
        for chunk = 1:length(Meno_chunks)
            [rho, theta] = CalculateAverageHeading(ftT,minVel, Meno_chunks{chunk,1});
            if ~isnan(theta)
                Meno_sum{count,1} = Meno_chunks{chunk,1}; 
                Meno_sum{count,2} = theta; 
                Meno_sum{count,3} = rho;
                count = count + 1; 
            end
        end
        count = 1; 
        for chunk = 1:length(not_Meno_chunks)
            [rho, theta] = CalculateAverageHeading(ftT,minVel, not_Meno_chunks{chunk,1});
            if ~isnan(theta)
                not_Meno_sum{count,1} = not_Meno_chunks{chunk,1}; 
                not_Meno_sum{count,2} = theta; 
                not_Meno_sum{count,3} = rho;
                count = count + 1; 
            end
        end
    

        % collect chunks with sim goals
        while ~isempty(Meno_sum)
            g = Meno_sum{1,2}; 
            other_g = [Meno_sum{2:end,2}];
            diffg = abs(wrapToPi(g - other_g));
            simg = find(diffg < pi/6) + 1; % allow a +/- 30 degree diff around 
            MenoData.Indices{mCount} = Meno_sum{1,1};%[Meno_sum{[1,simg],1}];
            [rho, theta] = CalculateAverageHeading(ftT,minVel, MenoData.Indices{mCount});
            speed = sqrt(ftT.velSide{1}(MenoData.Indices{mCount}).^2 + ftT.velFor{1}(MenoData.Indices{mCount}).^2);
            timeMov = sum(speed > 1.5)/60;
            MenoData.Goal(mCount) = theta; 
            MenoData.rho(mCount) = rho;
            MenoData.timeMov(mCount) = timeMov; 
            MenoData.Folder(mCount) = folder; 
            MenoData.numTrial(mCount) = nTrial;
            Meno_sum(1,:) = []; %Meno_sum([1, simg],:) = [];
            mCount = mCount + 1;
        end

        while ~isempty(not_Meno_sum)
            g = not_Meno_sum{1,2}; 
            other_g = [not_Meno_sum{2:end,2}];
            diffg = abs(wrapToPi(g - other_g));
            simg = find(diffg < pi/6) + 1; % allow a +/- 30 degree diff around 
            nMenoData.Indices{nCount} = not_Meno_sum{1,1};%[not_Meno_sum{[1,simg],1}];
            [rho, theta] = CalculateAverageHeading(ftT,minVel, nMenoData.Indices{nCount});
            speed = sqrt(ftT.velSide{1}(nMenoData.Indices{nCount}).^2 + ftT.velFor{1}(nMenoData.Indices{nCount}).^2);
            timeMov = sum(speed > 1.5)/60;
            nMenoData.Goal(nCount) = theta; 
            nMenoData.rho(nCount) = rho;
            nMenoData.timeMov(nCount) = timeMov;
            nMenoData.Folder(nCount) = folder;
            nMenoData.numTrial(nCount) = nTrial;
            not_Meno_sum(1,:) = []; %not_Meno_sum([1, simg],:) = [];
            nCount = nCount + 1;
        end
        
    end

    %% 


    end
end
