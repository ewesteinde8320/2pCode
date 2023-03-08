folder = 'Z:\2photon_data\PFL2_3_processedData\20220411\20220411-2_PFL2_fly1_FB'; 
nTrial = 1;

expID = get_expID(folder);
expList = {expID};

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
    load(fullfile(processedData_dir,['zscored_df_f_Trial_00',num2str(nTrial),'.mat']));
    load(fullfile(processedData_dir,['df_f_Trial_00',num2str(nTrial),'.mat']));
    load(fullfile(processedData_dir,['bump_parameters_Trial00',num2str(nTrial),'.mat']));
end


dff = [];
for roi = 1:size(dffData,1)
    dff(roi,:) = dffData.dff{roi};
end

Z = [];
for roi = 1:size(ZData,1)
    Z(roi,:) = ZData.Z{roi};
end

ftT.trialTime{1} = seconds(ftT.trialTime{1}); 
%%
timeStart = 320;
timeEnd = 420; 



vf = smoothdata(ftT.velFor{1}((ftT.trialTime{1} < timeEnd) & (ftT.trialTime{1} > timeStart)),'gaussian',50);
vy = smoothdata(ftT.velYaw{1}((ftT.trialTime{1} < timeEnd) & (ftT.trialTime{1} > timeStart)),'gaussian',50);
adj_rs = bump_params.adj_rs((ftT.trialTime{1} < timeEnd) & (ftT.trialTime{1} > timeStart));
pos_rad = bump_params.pos_rad((ftT.trialTime{1} < timeEnd) & (ftT.trialTime{1} > timeStart));
pos_rad = pos_rad(adj_rs > 0);
amp = bump_params.amp((ftT.trialTime{1} < timeEnd) & (ftT.trialTime{1} > timeStart));
amp = pos_rad(adj_rs > 0);
angle = ftT.cueAngle{1}((ftT.trialTime{1} < timeEnd) & (ftT.trialTime{1} > timeStart));
time = ftT.trialTime{1}((ftT.trialTime{1} < timeEnd) & (ftT.trialTime{1} > timeStart));

dff_chunk = dff(:,(ftT.trialTime{1} < timeEnd) & (ftT.trialTime{1} > timeStart));
Z_chunk = Z(:,(ftT.trialTime{1} < timeEnd) & (ftT.trialTime{1} > timeStart));
%%

fig1 = figure('Renderer', 'painters', 'Position', [100 300 1500 700],'color','w');
h1=subplot(2,1,1);
imagesc(dff_chunk)
colormap(flipud(gray))
title('dff activity');
set(gca,'YDir','normal')
box off
h2=subplot(2,1,2);
imagesc(Z_chunk)
colormap(flipud(gray))
title('Z activity');
set(gca,'YDir','normal')
box off

%%

fig2 = figure('Renderer', 'painters', 'Position', [100 300 1500 700],'color','w');
h1=subplot(5,1,1);
imagesc(Z_chunk)
colormap(flipud(gray))
colorbar
set(gca,'YDir','reverse')
h1.XAxis.Visible = 'off';
h1.YAxis.Visible = 'off';
box off
h2=subplot(5,1,2);
a = plot(time(adj_rs > 0),-wrapTo180(rad2deg(pos_rad + pi)),'k');
a.YData(abs(diff(a.YData)) > 180) = nan;
ylabel('bump pos')
xlim([min(time),max(time)])
ylim([-180,180])
set(gca, 'XTickLabel', [])
box off
h3=subplot(5,1,4);
b = plot(time,-angle,'k');
b.YData(abs(diff(b.YData)) > 200) = nan;
ylabel('head direction')
xlim([min(time),max(time)])
ylim([-180,180])
set(gca, 'XTickLabel', [])
box off
h4=subplot(5,1,3);
plot(time,amp,'k')
ylabel('bump amplitude')
xlim([min(time),max(time)])
set(gca, 'XTickLabel', [])
box off
h5=subplot(5,1,5);
plot(time,vf,'k')
ylabel('forward vel')
xlabel('Time (s)')
xlim([min(time),max(time)])
box off
linkaxes([h2,h3,h4,h5],'x')

%%
fig3 = figure('Renderer', 'painters', 'Position', [100 300 1500 700],'color','w');
h1=subplot(4,1,1);
imagesc(dff_chunk)
colormap(flipud(gray))
title('dff activity');
set(gca,'YDir','normal')
box off
h2=subplot(4,1,2);
plot(time,-angle,'k')
ylabel('head direction')
xlim([min(time),max(time)])
box off
h3=subplot(4,1,3);
plot(time,vy,'k')
ylabel('rotational vel')
xlim([min(time),max(time)])
box off
h4=subplot(4,1,4);
plot(time,vf,'k')
ylabel('forward vel')
xlim([min(time),max(time)])
box off
linkaxes([h2,h3,h4],'x')