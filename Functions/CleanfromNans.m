function spikes = CleanfromNans(spikes)
% Remove NaNs from events in dataset
k=1;
if ~isempty(spikes)
    while k <= size(spikes,1)
        if isnan(spikes(k,1)) || isnan(spikes(k,2))
            spikes(k,:) = [];
        else
            k = k+1;
        end
    end
end
end

function spikes = CleanOutliers(spikes)
% Create a rough polygon around AFG from which to exclude spikes (this
% function may be replaced with inpolygon for more accurate removal)
k=1;
if ~isempty(spikes)
    while k <= size(spikes,1)
        if spikes(k,1) < 0 || spikes(k,1) > 30 || spikes(k,2)<0 || spikes(k,2) > 30 ...
                || spikes(k,2) - 1.75*spikes(k,1) < -30
            spikes(k,:) = [];
        else
            k = k+1;
        end
    end
end
end
