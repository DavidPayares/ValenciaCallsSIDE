function spikes = CleanOutliers(spikes)
% Create a rough polygon around AFG from which to exclude spikes (this
% function may be replaced with inpolygon for more accurate removal)
k=1;
if ~isempty(spikes)
    while k <= size(spikes,1)
        if spikes(k,1) < 0 || spikes(k,1) > 36 || spikes(k,2)<0 || spikes(k,2) > 36
            spikes(k,:) = [];
        else
            k = k+1;
        end
    end
end
end
