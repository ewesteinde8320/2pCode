function [rho, theta, rho_all, theta_all] = CalculateAverageHeading_wholeSegment(ftT,window,minVel,data_idx, sampRate)

if strcmp(data_idx,'all')
        data_idx = 1:length(ftT.cueAngle{1});
end


behaviourData.angle = ftT.cueAngle{1}(data_idx); 
behaviourData.vel_side = ftT.velSide{1}(data_idx); 
behaviourData.vel_for = ftT.velFor{1}(data_idx);
behaviourData.vel_yaw = ftT.velYaw{1}(data_idx);
count = 1; 
mean_headingVectors = [];
idx_windows = [];
speed = abs(ftT.velFor{1}) + abs(ftT.velSide{1}) + abs(ftT.velYaw{1}*4.5);
window = window * sampRate; 
for i = 1 - window/2:1:length(behaviourData.angle) - window/2 
    idx = i:i + window; 
    if idx(end) > length(behaviourData.angle)
        idx = idx(1):1:length(behaviourData.angle);
    elseif idx(1) < 1 
        idx = 1:idx(end); 
    end
    angle_temp = behaviourData.angle(idx); 
    speed_temp = speed(idx); 
    angles_flyFor = angle_temp(speed_temp > minVel); 
    if ~isempty(angles_flyFor)
        x = cosd(angles_flyFor); 
        y = sind(angles_flyFor); %my arena has - angles to the left of the fly, + to the right, multiply y component by -1 to align physical arena coords to polar plot angles
        idx_windows(count,1) = idx(1);
        idx_windows(count,2) = idx(end);
        idx_windows(count,3) = length(find(speed_temp > minVel))/length(speed_temp)*100;
        mean_headingVectors(1,count)= sum(x)/length(x); 
        mean_headingVectors(2,count)= sum(y)/length(y); 
        count = count + 1; 
    else
        mean_headingVectors(1,count)= nan; 
        mean_headingVectors(2,count)= nan; 
        count = count + 1; 
    end
    
end 

trial_angle = behaviourData.angle(speed > minVel);
xTrial = cosd(trial_angle);
yTrial = sind(trial_angle); 
trialHeadingVector(1) = sum(xTrial)/length(xTrial);
trialHeadingVector(2) = sum(yTrial)/length(yTrial);

rho = sqrt(mean_headingVectors(1,:).^2 + mean_headingVectors(2,:).^2); 
theta = atan2(mean_headingVectors(2,:),mean_headingVectors(1,:)); 

rho_all = sqrt(trialHeadingVector(1).^2 + trialHeadingVector(2).^2); 
theta_all = atan2(trialHeadingVector(2),trialHeadingVector(1));