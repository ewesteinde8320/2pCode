function activityVSbehaviour_lineplots(ftT_down, Z, Zf, roiData, nTrial, expMd,savePlots,lineplotDir, expID)

    sum_mean = cell(3,1); 
    vy = wrapTo180((ftT_down.yawSpeed{1}/ (2*pi) ) * 360); 
    edges_vf = [min(ftT_down.fwSpeed{1}):0.5:max(ftT_down.fwSpeed{1})];
    edges_vs = [min(ftT_down.sideSpeed{1}):0.5:max(ftT_down.sideSpeed{1})];
    edges_vy = [min(vy):15:max(vy)];
    edges_angle = [-180:10:180];
    angle_down = ftT_down.cueAngle{1}; 
    
    trial_roiData = roiData(roiData.trialNum == nTrial,:);

    for run = 1:2
        sum_mean{1} = zeros(length(edges_vf)-1,1);
        sum_mean{2} = zeros(length(edges_vs)-1,1);
        sum_mean{3} = zeros(length(edges_vy)-1,1);
        sum_mean{4} = zeros(length(edges_angle)-1,1);

        if run == 1
            activityTable = Z;
            figure(Name=['Zscore vs behaviour, trial ', num2str(nTrial)]);clf
            label = 'Z'; 
        else
            activityTable = Zf;
            figure(Name=['zf_f vs behaviour, trial ', num2str(nTrial)]);clf
            label = 'df_f';
        end

        if ~contains(expMd.expName{1},'LAL')
            for roi = 1:size(trial_roiData,1)
                activity = activityTable.data(activityTable.roiName == roi);

                % vf

                behaviour = ftT_down.fwSpeed{1}; 
                [vf_zscore, centers_vf] = binData(activity, behaviour, edges_vf);
                sum_mean{1} = sum_mean{1} + vf_zscore; 

                subplot(4,1,1);
                plot(centers_vf,vf_zscore)
                colororder(parula(roi))
                ylabel(label)
                xlabel('vf (mm/s)')
                hold on


                % vs 
                behaviour = ftT_down.sideSpeed{1}; 
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
                [angle_zscore, centers_angle] = binData(activity, angle_down, edges_angle);
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
                    % vf
                    behaviour = ftT_down.fwSpeed{1}; 
                    [vf_zscore, centers_vf] = binData(activity, behaviour, edges_vf);
                    sum_mean{1}(:,roi) = vf_zscore; 

                    l(1) = subplot(4,2,1);
                    plot(centers_vf,vf_zscore)
                    ylabel(label)
                    xlabel('vf (mm/s)')
                    legend(trial_roiData.roiName,'Interpreter', 'none');
                    hold on


                    % vs 
                    behaviour = ftT_down.sideSpeed{1}; 
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
                    [angle_zscore, centers_angle] = binData(activity, ftT_down.cueAngle{1}, edges_angle);
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
        saveas(z, fullfile(lineplotDir,[expID,'_',num2str(nTrial),'_zScore_behaviour_lineplots.fig']));
        saveas(zf, fullfile(lineplotDir,[expID,'_',num2str(nTrial),'_df_f_behaviour_lineplots.fig']));
    end
end