%%
% 
% Find the fs lower-bound considering CP length and subframe duration.
% Return NaN if the lower-bound is invalid.
% 
function fs_LB = find_fs_LB(gfdm, Channel, fs_UB)

% CP length for fs=Channel
CP = ceil(gfdm.rms*Channel*1e-3);

if Channel >= (gfdm.D + CP)/gfdm.tau_subframe
    % If fs=Channel fits the subframe limit
    fs_LB = Channel;
elseif gfdm.rms == 0
    % If fs=Channel DOESN'T fit the subframe limit, and RMS_delay=0
    fs_LB = gfdm.D/gfdm.tau_subframe;
    if fs_LB > fs_UB
        fs_LB = nan;
    end
else
    % If fs=Channel DOESN'T fit the subframe limit, and RMS_delay!=0
    % Need to find fs_LB

    % Find fs_LB to the first digit
    fs_LB = Channel;
    while 1
        fs_right = (gfdm.D + ceil(gfdm.rms*fs_LB*1e-3))/gfdm.tau_subframe;
        if fs_LB > fs_UB
            fs_LB = nan;
            return;
        elseif fs_LB >= fs_right
            break;
        end
        fs_LB = fs_LB+1;
    end

    if fs_LB == fs_right
        return;
    end
    
    % Find fs_LB to the first decimal
    fs_LB = fs_LB-1;
    while 1
        fs_LB = fs_LB + 0.1;
        fs_right = (gfdm.D + ceil(gfdm.rms*fs_LB*1e-3))/gfdm.tau_subframe;
        if fs_LB >= fs_right
            break;
        end
    end

    if fs_LB == fs_right
        return;
    end

    % Find fs_LB to the second decimal
    fs_LB = fs_LB-0.1;
    while 1
        fs_LB = fs_LB + 0.01;
        fs_right = (gfdm.D + ceil(gfdm.rms*fs_LB*1e-3))/gfdm.tau_subframe;
        if fs_LB >= fs_right
            break;
        end
    end
end

end