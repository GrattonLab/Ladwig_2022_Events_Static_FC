clear all;
% sample data in 4 diff ways and see the effects on
% network structure like in Fig 3

sessions = 1:10;
subjects = {'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07','MSC09', 'MSC10'};

%% 

for k = 1:length(subjects)
    disp(k)
    subject = subjects{k};
    for j = 1:length(sessions)
        data_path = 'YOUR DATA PATH';
        thresh = 0.05;
        bins = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, thresh];
        bin = bins(1);

        session = j;
        functionals_path = [data_path subject '_parcel_timecourse.mat'];
        load(functionals_path)
        parcel_timeseries = parcel_time{session};
        tmask = tmask_all{session};
        timeseries = parcel_timeseries(tmask==1, :);
        fullfc = corr(timeseries);
        numpts = round(thresh*size(timeseries,1));
        startpoint = floor((1-bin)*size(timeseries,1))+1;
        endpoint = startpoint + numpts-1;

        % for each session we sample randomly, high, low, and consectuively
        [keepfc, events_idx, rms_idx, timeseries_z] = getEventsFc(timeseries, startpoint, endpoint);
        rand_idx = datasample(rms_idx, numpts, 'Replace', false);
        ts_random = timeseries_z(rand_idx,:);
        fc_random = corr(ts_random);

        ts_high = timeseries_z(rms_idx(1:numpts), :);
        fc_high = corr(ts_high);
        ts_low = timeseries_z(rms_idx(length(rms_idx)-numpts+1: length(rms_idx)), :);
        fc_low = corr(ts_low);

% as to not bias the results, we iterate over all the possible starting
% locations for consecutive comparisons and then later take the mean of
% them
        for i = 1:length(rms_idx)
            start_consec = i;
            end_consec = start_consec + numpts;
            if end_consec > length(rms_idx)
                pre = start_consec:length(rms_idx);
                post = 1:(end_consec - length(rms_idx)-1);
                idx = [pre post];
                ts_random_consec = timeseries_z(idx, :);

            else 
                ts_random_consec = timeseries_z(start_consec:start_consec+numpts-1, :);
            end
            newfc = corr(ts_random_consec);
            compare_fc(i) = corr(newfc(:), fullfc(:));
        end

        compare(k, j, 1) = corr(fc_random(:), fullfc(:));
        compare(k, j, 2) = corr(fc_low(:), fullfc(:));
        compare(k, j, 3) = corr(fc_high(:), fullfc(:));
        compare(k, j, 4) = mean(compare_fc);
    end
end


%% function 
function [keepfc, events_idx, rms_idx, timeseries_z] = getEventsFc(timeseries, startpoint, endpoint)
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

    eventCofluxRms = sqrt(sum(coflux.^2));
    [~, rms_idx] = sort(eventCofluxRms, 'descend');
    events_idx = rms_idx(startpoint:endpoint);
    keepfc = corr(timeseries_z(rms_idx(startpoint:endpoint),:));
end