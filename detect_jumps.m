function jump_array = detect_jumps(roiData, ftT_down, jump_window)

    yChannel = ftT_down.PanelsY{1}; 
    yChannel_temp = yChannel - [1 2 3 4];
    
    for row = 1:size(yChannel_temp,1)
        yChannel_idx = knnsearch(yChannel_temp(row,:)',0);
        yChannel_clean(row,1) = yChannel_idx;
    end
    
    yChannel = yChannel_clean; 

    count = 1; 
    prev_jump = 0; 
    transition = zeros(size(yChannel)); 

    for idx = 1:length(yChannel)
        value = yChannel(idx);
        if isnan(value)
            print ('nan approaching')
        end

        if count ~= 1 
            if value ~= prev_val && ~(isnan(value)) && (count - prev_jump > 5) 
                transition(count) = 1;
                prev_jump = count; 
            else
                transition(count) = 0; 
            end

            if ~(isnan(value))
                prev_val = value;
            else
                prev_val = prev_val;
            end
            
            count = count + 1; 

        else
            transition(count) = 0; 
            prev_val = value;
            count = count + 1;
        end
    end

    volRate = roiData.sampRate(1); 
    window_datapoints = round(jump_window*volRate); 
    jump_idx = find(transition == 1); 
    jump_array = zeros(length(jump_idx),3); 

    count = 1; 

    for idx = 1:length(jump_idx)
        pre_jump = jump_idx(idx) - window_datapoints;
        post_jump = jump_idx(idx) + window_datapoints;
        jump_array(count,1) = pre_jump; 
        jump_array(count,2) = jump_idx(idx); 
        jump_array(count,3) = post_jump; 
        count = count + 1;
    end
    
    %Debugging plot 
%     figure()
%     plot(yChannel)
%     hold on
%     plot(transition,'ro')
    
end