%% Code to calculate cofluctuation and bin timepoints for all MSC subjects as in Fig 1
clear all;
subjects = {'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07','MSC09', 'MSC10'};
sessions = 1:10;
session_type = 'rest';
thresh = 0.05;
bins = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, thresh];
mod_threshold = 0.05;

%% this is good for the coflux versus network structure measures 
for i = 1:length(bins)
    disp(i) 
    for j = 1:length(subjects)
        for k = 1:length(sessions)
            functionals_path = ['YOUR FILE PATH HERE' session_type '/' subjects{j} '_parcel_timecourse.mat'];
            load(functionals_path)
            parcel_timeseries = parcel_time{sessions(k)};
            tmask = tmask_all{sessions(k)};
            timeseries = parcel_timeseries(tmask==1, :);

            % filter out sessions with < 333 data pts 
            if(size(timeseries,1) <= 333)
                fc(i,j,k) = NaN;
                mod(i,j,k) = NaN;
                continue;
            end

            fullfc = corr(timeseries);
           
            numpts = round(thresh*size(timeseries,1));
            startpoint = floor((1-bins(i))*size(timeseries,1))+1;
            endpoint = startpoint + numpts -1;

            % for every subject/session/bin get an FC matrix and modularity
            % value
            keepfc = getEventsFc(timeseries, startpoint, endpoint);
            fc(i, j, k) = corr(fullfc(:), keepfc(:));
            [mat_thresh r kden] = matrix_thresholder(keepfc, mod_threshold ,'kden');
            mat_thresh(find(mat_thresh > 0))=1;
            [Ci Q] = modularity_und(mat_thresh);
            mod(i, j, k) = Q;
        end 
    end
end

%% you may want to save these out to use them since they take a while to generate
save('discreteness_output.mat', 'fc', 'mod')

%% function 
function [keepfc, events_idx, timeseries_z] = getEventsFc(timeseries, startpoint, endpoint)
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