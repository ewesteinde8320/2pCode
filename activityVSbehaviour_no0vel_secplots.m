function activityVSbehaviour_no0vel_secplots(ftT_down, Z, Zf, roiData, nTrial, threshold, expMd,savePlots,lineplotDir, expID)

%% Remove idx where the fly isn't moving
    
    %total_mov_mm = abs(fft_data[:,vel_for_idx]) + abs(fft_data[:,vel_side_idx]) + abs(fft_data[:,vel_yaw_idx]*4.5)
    no0vel_idx = find(ftT_down.moveSpeed{1} > threshold);
    %speed = sqrt(ftT_down.sideSpeed{1}.^2 + ftT_down.fwSpeed{1}.^2);
    %vy = ftT_down.yawSpeed{1}; 
    %no0vel_idx = find(abs(ftT_down.sideSpeed{1}) > 0.5);
    %no0vel_idx = find(speed > threshold); 
    vf = ftT_down.fwSpeed{1}(no0vel_idx);
    vs = ftT_down.sideSpeed{1}(no0vel_idx);
    vy = ftT_down.yawSpeed{1}(no0vel_idx);
    angle = ftT_down.cueAngle{1}(no0vel_idx); 
    
%%

    sum_mean = cell(3,1); 
    vy = wrapTo180((vy/ (2*pi) ) * 360); 
    edges_vf = [min(vf):0.5:max(vf)];
    edges_vs = [min(vs):0.5:max(vs)];
    edges_vy = [min(vy):10:max(vy)];
    edges_angle = [-180:10:180]; 
    
    trial_roiData = roiData(roiData.trialNum == nTrial,:);

    for run = 1:2
        sum_mean{1} = zeros(length(edges_vf)-1,1);
        sum_mean{2} = zeros(length(edges_vs)-1,1);
        sum_mean{3} = zeros(length(edges_vy)-1,1);
        sum_mean{4} = zeros(length(edges_angle)-1,1);

        if run == 1
            activityTable = Z;
            figure(Name=['Zscore vs behaviour no 0 vel , trial ', num2str(nTrial)]);clf
            label = 'Z'; 
        else
            activityTable = Zf;
            figure(Name=['zf_f vs behaviour, trial no 0 vel', num2str(nTrial)]);clf
            label = 'df_f';
        end

        if ~contains(expMd.expName{1},'LAL')
            for roi = 1:size(trial_roiData,1)
                activity = activityTable.data(activityTable.roiName == roi);
                activity = activity(no0vel_idx);

                % vf

                behaviour = vf; 
                [vf_zscore, centers_vf] = binData(activity, behaviour, edges_vf);
                sum_mean{1} = sum_mean{1} + vf_zscore; 

                subplot(4,1,1);
                plot(centers_vf,vf_zscore)
                colororder(parula(roi))
                ylabel(label)
                xlabel('vf (mm/s)')
                hold on


                % vs 
                behaviour = vs; 
                [vs_zscore, centers_vs] = binData(activity, behaviour, edges_vs);
                sum_mean{2} = sum_mean{2} + vs_zscore; 

                subplot(4,1,2);
                plot(centers_vs,vs_zscore)
                colororder(parula(roi))
                ylabel(label)
                xlabel('vs (mm/s)')
                hold on

                % vy 
                [vy_zscore, centers_vy] = binData(activity, vy, edges_vy);
                sum_mean{3} = sum_mean{3} + vy_zscore; 


                subplot(4,1,3);
                plot(centers_vy,vy_zscore)
                colororder(parula(roi))
                ylabel(label)
                xlabel('vy (deg/s)')
                hold on

                % angle
                [angle_zscore, centers_angle] = binData(activity, angle, edges_angle);
                sum_mean{4} = sum_mean{4} + angle_zscore; 


                subplot(4,1,4);
                plot(centers_angle,angle_zscore)
                colororder(parula(roi))
                ylabel(label)
                xlabel('cue pos (deg)')
                hold on
            end


                subplot(4,1,1)
                plot(centers_vf,sum_mean{1}/size(trial_roiData,1),'k','LineWidth',1.5)
                subplot(4,1,2)
                plot(centers_vs,sum_mean{2}/size(trial_roiData,1),'k','LineWidth',1.5)
                subplot(4,1,3)
                plot(centers_vy,sum_mean{3}/size(trial_roiData,1),'k','LineWidth',1.5)
                subplot(4,1,4)
                plot(centers_angle,sum_mean{4}/size(trial_roiData,1),'k','LineWidth',1.5)

        else
               for roi = 1:size(trial_roiData,1)
                    activity = activityTable.data(activityTable.roiName == roi);
                    activity = activity(no0vel_idx);
                    
                    % vf
                    behaviour = vf; 
                    [vf_zscore, centers_vf] = binData(activity, behaviour, edges_vf);
                    sum_mean{1}(:,roi) = vf_zscore; 

                    l(1) = subplot(4,2,1);
                    plot(centers_vf,vf_zscore)
                    ylabel(label)
                    xlabel('vf (mm/s)')
                    legend(trial_roiData.roiName,'Interpreter', 'none');
                    hold on


                    % vs 
                    behaviour = vs; 
                    [vs_zscore, centers_vs] = binData(activity, behaviour, edges_vs);
                    sum_mean{2}(:,roi) = vs_zscore; 

                    l(2) = subplot(4,2,3);
                    plot(centers_vs,vs_zscore)
                    ylabel(label)
                    xlabel('vs (mm/s)')
                    hold on

                    % vy 
                    [vy_zscore, centers_vy] = binData(activity, vy, edges_vy);
                    sum_mean{3}(:,roi) = vy_zscore; 


                    l(3) = subplot(4,2,5);
                    plot(centers_vy,vy_zscore)
                    ylabel(label)
                    xlabel('vy (deg/s)')
                    hold on

                    % angle
                    [angle_zscore, centers_angle] = binData(activity, angle, edges_angle);
                    sum_mean{4}(:,roi) = angle_zscore; 


                    l(4) = subplot(4,2,7);
                    plot(centers_angle,angle_zscore)
                    ylabel(label)
                    xlabel('cue pos (deg)')
                    hold on 
               end

                subplot(4,2,2)
                plot(centers_vf,sum_mean{1}(:,1)-sum_mean{1}(:,2))
                title('L-R')
                
                subplot(4,2,4)
                plot(centers_vs,sum_mean{2}(:,1)-sum_mean{2}(:,2))

                subplot(4,2,6)
                plot(centers_vy,sum_mean{3}(:,1)-sum_mean{3}(:,2))

                subplot(4,2,8)
                plot(centers_angle,sum_mean{4}(:,1)-sum_mean{4}(:,2))

        end

        if run == 1
            z = gcf;
        else
            zf = gcf;
        end
    end


    if savePlots == 1
        saveas(z, fullfile(lineplotDir,[expID,'_',num2str(nTrial),'_zScore_behaviour_no0vel.fig']));
        saveas(zf, fullfile(lineplotDir,[expID,'_',num2str(nTrial),'_df_f_behaviour_no0vel.fig']));
    end
end