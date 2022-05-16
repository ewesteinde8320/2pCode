function plot_jumps(jump_array, jump_array_down, Z, ftData, expID, nTrial, jumpDir, roiData)

%% plot jump windows


    figure()
    for jump = 1:length(jump_array)
        window_idx = [jump_array(jump,1):jump_array(jump,3)];
        jump_cueAngle = ftData.cueAngle{1}(window_idx);
        jump_idx = round(length(jump_cueAngle)/2);
        time = [-windowSize:1/60:windowSize];

        subplot(length(jump_array),1,jump)
        plot(time,jump_cueAngle)
        %xline(0)
        ylim([-200 200])
        title('Cue Angles')
    end

    if savePlots == 1
        saveas(gcf, fullfile(jumpDir,[expID,'_',num2str(nTrial),'_cueAngleSum.fig']));
    end


    figure()
    for jump = 1:length(jump_array_down)
        window_idx = [jump_array_down(jump,1):jump_array_down(jump,3)];
        time = [-windowSize:1/roiData.sampRate(1):windowSize];
        subplot(length(jump_array_down),1,jump)
        for roi = 1:length(unique(Z.roiName))
            jump_Z = Z(Z.roiName == roi,:).data(window_idx);
            roi_jump(roi,:) = jump_Z; 
            jump_idx = round(length(jump_Z)/2);
            plot(time,jump_Z)
            hold on
        end
        jump_aveZ = sum(roi_jump,1)/length(unique(Z.roiName));
        plot(time,jump_aveZ,'k','lineWidth',2)
        xline(0)
        title('Z')
    end

        if savePlots == 1
            saveas(gcf, fullfile(jumpDir,[expID,'_',num2str(nTrial),'_roiDataSum.fig']));
        end

    %% individual jump plots


    for jump = 1:length(jump_array)
        i = figure();
        window_idx = [jump_array(jump,1):jump_array(jump,3)];
        jump_cueAngle = ftData.cueAngle{1}(window_idx);
        jump_idx = round(length(jump_cueAngle)/2);
        time = [-windowSize:1/60:windowSize];

        subplot(14,1,1:4)
        plot(time,jump_cueAngle)
        %xline(0)
        ylim([-200 200])
        ylabel('cue angle')
        title(['jump ', num2str(jump)])

        subplot(14,1,9:10)
        plot(time, smoothdata(ftData.fwSpeed{1}(window_idx),'loess',20))
        ylabel('vf')
        xline(0)
        subplot(14,1,11:12)
        plot(time, smoothdata(ftData.yawSpeed{1}(window_idx),'loess',20))
        ylabel('vy')
        xline(0)
        subplot(14,1,13:14)
        plot(time, smoothdata(ftData.sideSpeed{1}(window_idx),'loess',20))
        ylabel('vs')
        xline(0)


        window_idx = [jump_array_down(jump,1):jump_array_down(jump,3)];
        time = [-windowSize:1/roiData.sampRate(1):windowSize];
        subplot(14,1,5:8)
        for roi = 1:length(unique(Z.roiName))
            jump_Z = Z(Z.roiName == roi,:).data(window_idx);
            roi_jump(roi,:) = jump_Z; 
            jump_idx = round(length(jump_Z)/2);
            plot(time,jump_Z)
            hold on
        end
        jump_aveZ = sum(roi_jump,1)/length(unique(Z.roiName));
        plot(time,jump_aveZ,'k','lineWidth',2)
        xline(0)
        ylabel('Z')

    %     subplot(14,1,9:10)
    %     plot(time, ftT_down.fwSpeed{1}(window_idx))
    %     ylabel('vf')
    %     xline(0)
    %     subplot(14,1,11:12)
    %     plot(time, ftT_down.yawSpeed{1}(window_idx))
    %     ylabel('vy')
    %     xline(0)
    %     subplot(14,1,13:14)
    %     plot(time, ftT_down.sideSpeed{1}(window_idx))
    %     ylabel('vs')
    %     xline(0)

        if savePlots == 1
            saveas(i, fullfile(jumpDir,[expID,'_',num2str(nTrial),'_jump',num2str(jump),'.fig']));
        end

    end

end
    