function [triggerIdx, rho, Meno_chunks, not_Meno_chunks] = SegmentMenovsNotMeno(ftT, window, minVel,highThres,lowThres, folder, trial,toPlot)

triggerIdx = [];
rho = [];
Meno_chunks = {}; 
not_Meno_chunks = {};

try
    
        sampRate = 60;     
        window = window * sampRate; 
        if window > length(ftT.velFor{1})
            window = length(ftT.velFor{1}); 
        end
        count = 1; 


        mean_headingVectors = [];
        speed = sqrt(ftT.velFor{1}.^2 + ftT.velSide{1}.^2);
        
        [~, jump_array, ~] = detect_jumps(ftT, 5, 5,0);
        jump_idx = [];
        for jump = 1:size(jump_array,1)
            jump_idx = [jump_idx , jump_array(jump,2):jump_array(jump,3)];
        end

        for i = 1 - window/2:1:length(ftT.cueAngle{1}) - window/2 
            idx = i:i + window; 
            if idx(end) > length(ftT.cueAngle{1})
                idx = idx(1):1:length(ftT.cueAngle{1});
            elseif idx(1) < 1 
                idx = 1:idx(end); 
            end

            if sum(ismember(idx, jump_idx)) % remove idx from window that were influenced by jumps
                badIdx = ismember(idx, jump_idx);
                idx = idx(badIdx == 0);
            end

                angle_temp = ftT.cueAngle{1}(idx); 
                speed_temp = speed(idx); 
                angles_flyFor = angle_temp(speed_temp > minVel); 
                if ~isempty(angles_flyFor) 
                    x = cosd(angles_flyFor); 
                    y = sind(angles_flyFor); 
                    idx_windows{count,1} = idx;
                    mean_headingVectors(1,count)= sum(x)/length(x); 
                    mean_headingVectors(2,count)= sum(y)/length(y);
                    count = count + 1;
                else
                    mean_headingVectors(1,count)= nan; 
                    mean_headingVectors(2,count)= nan; 
                    slope(count) = nan;
                    count = count + 1;
                end
        end

        rho = sqrt(mean_headingVectors(1,:).^2 + mean_headingVectors(2,:).^2); 
        
    %% sinuosity index
    
%     [xPos, yPos] = plot2DTrajectory(folder, 0, trial);
%     count = 1; 
% 
%  for i = 1 - window/2:1:length(ftT.cueAngle{1}) - window/2
%             idx = i:i + window; 
%             if idx(end) > length(ftT.cueAngle{1})
%                 idx = idx(1):1:length(ftT.cueAngle{1});
%             elseif idx(1) < 1 
%                 idx = 1:idx(end); 
%             end
%             
%             if sum(ismember(idx, jump_idx)) % remove idx from window that were influenced by jumps
%                 badIdx = ismember(idx, jump_idx);
%                 idx = idx(badIdx == 0);
%             end
%             
%             xChunk = xPos(idx);
%             yChunk = yPos(idx);
%             
%             xdiff = abs(gradient(xChunk)); 
%             ydiff = abs(gradient(yChunk));
%             
%             edist = sqrt((xChunk(end)-xChunk(1)).^2 + (yChunk(end)-yChunk(1)).^2);
%             curveLength = sum(sqrt(xdiff.^2 + ydiff.^2));
%             
%             ratio = curveLength/edist;
%             
% %             if ratio > 10
% %                 disp(ratio)
% %                 figure();patch('XData',xChunk,'YData',yChunk,'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
% %                 hold on
% %                 patch('XData',xChunk([1,end]),'YData',yChunk([1,end]),'EdgeColor','r','FaceColor','none','LineStyle','-','Marker','.', 'MarkerSize',7);
% %             pause
% %             end
% 
%             
%             sIdx(count) = -ratio;
%             
%             count = count + 1;     
%  end
 
%         figure();
%         ax1 = subplot(3,1,1);
%         plot(ftT.cueAngle{1})
%         ax2 = subplot(3,1,2);
%         plot(sIdx)
%         yline(-2)
%         yline(-3)
%         ax3 = subplot(3,1,3);
%         plot(rho)
%         yline(.88)
%         yline(0.8)
%         
%         
%        linkaxes([ax1,ax2,ax3],'x')

%% 

        triggerIdxRho = schmittTrigger(rho',highThres,lowThres);
        %triggerIdxSidx = schmittTrigger(sIdx',-2,-3);
        
        triggerIdx = zeros(size(triggerIdxRho)); 
        triggerIdx = triggerIdxRho;
        %triggerIdx(triggerIdxRho == 1 & triggerIdxSidx == 1) = 1;

        pastMenoIdx = 1;
        pastNotMenoIdx = 1;
        for i = 1:length(triggerIdx)
            iValue = triggerIdx(i);
            if i > 1
                if iValue == 1
                    timeFromMeno = i - pastMenoIdx; 
                    if timeFromMeno < 2 * sampRate % last transition from menotaxis occured less than 5 seconds before
                        triggerIdx(pastMenoIdx:i) = 1; 
                    end
                    pastMenoIdx = i; 
                elseif iValue == 0
                    timeFromNotMeno = i - pastNotMenoIdx;
                    if timeFromNotMeno < 2 * sampRate
                        triggerIdx(pastNotMenoIdx:i) = 0; 
                    end
                    pastNotMenoIdx = i; 
                end
            end
        end
        
        

        %% Group meno & not meno data chunks in arrays

            chunk = []; 
            Mcount = 1; 
            nMcount = 1; 
            for idx = 1:length(triggerIdx)
                if isempty(chunk)
                    if ~isnan(triggerIdx(idx))
                        chunk(1) = idx;
                        chunkValue = triggerIdx(idx);
                    end
                else
                    if triggerIdx(idx) == chunkValue || isnan(triggerIdx(idx))
                        chunk = [chunk, idx];
                    else
                        if chunkValue == 1 
                            Meno_chunks{Mcount} = chunk;
                            Mcount = Mcount + 1; 
                            chunk = [];
                        else
                            not_Meno_chunks{nMcount} = chunk;
                            nMcount = nMcount + 1; 
                            chunk = [];
                        end
                    end
                end
            end

     %% plot 2D path trajectory meno = red not meno = black
if toPlot

                    yawAngPos = ftT.cueAngle{1};
                    fwdAngVel = ftT.velFor{1};
                    slideAngVel = ftT.velSide{1};

                    sampRate = 60; % check
                    % conversion factor between degrees and mm
                    circum = 9 * pi; % circumference of ball, in mm
                    mmPerDeg = circum / 360; % mm per degree of ball

                    % position incorporating heading - as if fly were walking on x-y plane,
                    %  x-y coordinates at each time point
                    % start with fly at (0,0) and facing 0 deg
                    zeroedYawAngPos = yawAngPos - yawAngPos(1); 

                    % movement in x (in degrees) at each time point
                    xChangePos = (fwdAngVel ./ sampRate) .* sind(zeroedYawAngPos) + ...
                        (slideAngVel ./ sampRate) .* sind(zeroedYawAngPos + 90);  

                    % x position in mm (i.e. x-coordinate of fly's position at each time 
                    %  point), starts at 0
                    xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;
                    minX = min(xPos);
                    maxX = max(xPos);

                    % movement in y (in degrees) at each time point
                    yChangePos = (fwdAngVel ./ sampRate) .* cosd(zeroedYawAngPos) + ...
                        (slideAngVel ./ sampRate) .* cosd(zeroedYawAngPos + 90);

                    % y position in mm (i.e. y-coordinate of fly's position at each time 
                    %  point), starts at 0
                    yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;
                    minY = min(yPos);
                    maxY = max(yPos);

                    %[~, ~, transition] = detect_jumps(ftT, 2, 2);
                    time = seconds(ftT.trialTime{1});
                    nMenoTime = time(triggerIdx == 0);
                    MenoTime = time(triggerIdx ==1); 
                    timeStart = time(1); 
                    timeEnd = max(time); 
                    nxPos = xPos(triggerIdx == 0);
                    nyPos = yPos(triggerIdx == 0);
                    angle =  -ftT.cueAngle{1}; 
                    %angle(speed < minVel) = nan; 
                    nAngle = angle(triggerIdx ==0); 
                    %rho(speed < minVel) = nan;
                    nRho = rho(triggerIdx == 0); 
                    
                    mAngle = angle(triggerIdx ==1);
                    mxPos = xPos(triggerIdx == 1);
                    myPos = yPos(triggerIdx == 1);
                    mRho = rho(triggerIdx == 1); 
                    %jumpTimes = time(transition == 1); 
                    figure(66);
                    patch('XData',nxPos(nMenoTime > timeStart & nMenoTime < timeEnd),'YData',nyPos(nMenoTime > timeStart & nMenoTime < timeEnd),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
                    hold on
                    patch('XData',mxPos(MenoTime > timeStart &  MenoTime < timeEnd),'YData',myPos(MenoTime > timeStart &  MenoTime < timeEnd),'EdgeColor','r','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
                    %patch('XData',xPos(1),'YData',yPos(1),'EdgeColor','g','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',20);
                    %patch(xPos(transition == 1),yPos(transition == 1),time(transition == 1),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',7);
                    set(gcf,'color','w');
                    xlabel('mm')
                    xlim([min(minX,minY),max(maxX, maxY)])
                    ylim([min(minX,minY),max(maxX, maxY)])

                    %idxes = regexp(folder,'\');
                    title(['window: ',num2str(window/sampRate),' highThres: ', num2str(highThres),' lowThres: ', num2str(lowThres)], 'Interpreter', 'none')

                    saveas(gcf, fullfile(folder,'images',['menoVsnonMeno_Path_schmitt_sinuosity_00',num2str(trial),'.fig']))
                    
                    figure(77);clf;
                    set(gcf,'color','w','renderer','painters')
                    h(1) =  subplot(2,1,1);
                    hold on
                    a = plot(nMenoTime(nMenoTime > timeStart & nMenoTime < timeEnd),nAngle(nMenoTime > timeStart & nMenoTime < timeEnd),'k');
                    a.XData(abs(diff(a.XData)) > 10) = nan;
                    b = plot(MenoTime(MenoTime > timeStart & MenoTime < timeEnd),mAngle(MenoTime > timeStart & MenoTime < timeEnd),'r');
                    try
                        b.XData(abs(diff(b.XData)) > 10) = nan;
                    catch
                    end
                    ylabel('HD')
                    ylim([-180,180])
                    xlim([timeStart,timeEnd])
                    box off
                    h(2) = subplot(2,1,2);
                    hold on
                    c = plot(nMenoTime(nMenoTime > timeStart & nMenoTime < timeEnd),nRho(nMenoTime > timeStart & nMenoTime < timeEnd),'k');
                    c.XData(abs(diff(c.XData)) > 10) = nan;
                    d = plot(MenoTime(MenoTime > timeStart & MenoTime < timeEnd),mRho(MenoTime > timeStart & MenoTime < timeEnd),'r');
                    try
                        d.XData(abs(diff(d.XData)) > 10) = nan;
                    catch
                    end
                    %plot(time(time > timeStart & time < timeEnd),triggerIdx(time > timeStart & time < timeEnd))
                    ylabel('rho')
                    box off
                    xlim([timeStart,timeEnd])
                    linkaxes(h,'x')

    %     figure();
    %     ax1 = subplot(3,1,1);
    %     plot(seconds(ftT.trialTime{1}),ftT.cueAngle{1})
    %     xlim([min(seconds(ftT.trialTime{1})),max(seconds(ftT.trialTime{1}))])
    %     ax2 = subplot(3,1,2);
    %     plot(seconds(ftT.trialTime{1}),slope)
    %     yline(-30)
    %     yline(-100)
    %     ax3 = subplot(3,1,3);
    %     plot(seconds(ftT.trialTime{1}),triggerIdx)
    %     linkaxes([ax1, ax2, ax3],'x')
    %     
    %     figure();
    %     ax1 = subplot(4,1,1);
    %     plot(seconds(ftT.trialTime{1}),ftT.cueAngle{1})
    %     xlim([min(seconds(ftT.trialTime{1})),max(seconds(ftT.trialTime{1}))])
    %     ax2 = subplot(4,1,2);
    %     plot(seconds(ftT.trialTime{1}),slope)
    %     yline(0)
    %     yline(-100)
    %     xlim([min(seconds(ftT.trialTime{1})),max(seconds(ftT.trialTime{1}))])
    %     ax3 = subplot(4,1,3);
    %     plot(seconds(ftT.trialTime{1}),rho)
    %     xlim([min(seconds(ftT.trialTime{1})),max(seconds(ftT.trialTime{1}))])
    %     ax4 = subplot(4,1,4);
    %     plot(seconds(ftT.trialTime{1}),rad2deg(theta))
    %     xlim([min(seconds(ftT.trialTime{1})),max(seconds(ftT.trialTime{1}))])
    %     linkaxes([ax1, ax2, ax3, ax4],'x')
end
catch
    disp(['folder ',folder,'trial_',num2str(trial),' failed'])
end