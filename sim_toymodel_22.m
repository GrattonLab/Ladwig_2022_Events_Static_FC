clear all;
% goal of this script is to make the 2X2 matrix which serves as a toy model
% and do events analysis on it as in Fig 4 

%% make sample 
x = 0:pi/24:4*pi;
signal(1,:) = sin(x);

nodes = 2;
networks = 2;
node_names = {'A1', 'A2', 'B1', 'B2'};
for i = 1:nodes
    signal(i,:) = sin(x);
    signal(i+nodes,:) = sin(x+pi);
end

%% generate data a bunch of times 
for j = 1:1000
    noise = randn(size(signal));
    signal_noise(j, :, :) = signal + 0.5*noise;
    real_fc(j, :, :) = corr(squeeze(signal_noise(j,:,:))');

    bins = [1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.05];
    thresh = 0.05; 
    numpts = round(size(signal,2)*thresh);
 
    % do the events analysis on top of it 
    for i = 1:length(bins)
        startpoint = floor((1-bins(i))*size(signal,2))+1;
        endpoint = startpoint + numpts -1;
        [keepfc(j,i,:,:), eventCofluxRms(j,:)] = getEventsFc(squeeze(signal_noise(j,:,:))', startpoint, endpoint);
        clean_idx = find(tril(abs(squeeze(keepfc(j,i,:,:))),-1) > 0);
        clean_real_fc = real_fc(j,:,:);
        clean_keepfc = keepfc(j, i, :, :);
        fc(j, i) = corr(clean_real_fc(clean_idx), clean_keepfc(clean_idx));
    end
end

%% 
function [keepfc, eventCofluxRms] = getEventsFc(timeseries, startpoint, endpoint)

    timeseries_z = zscore(timeseries);
    [time,nodes] = size(timeseries_z);
    coflux = zeros(nodes*(nodes-1)/2,time);   
    count = 0;
    
    for i=1:nodes
        for j=i+1:nodes
            count = count+1;
            time1 = timeseries_z(:, i);
            time2 = timeseries_z(:, j);
            coflux(count,:) = time1.*time2;
        end
    end


    eventCofluxRms = sqrt(sum(coflux.^2, 1));
    [~, rms_idx] = sort(eventCofluxRms, 'descend');
    events_idx = rms_idx(startpoint:endpoint);
    keepfc = corr(timeseries_z(rms_idx(startpoint:endpoint),:));

end


