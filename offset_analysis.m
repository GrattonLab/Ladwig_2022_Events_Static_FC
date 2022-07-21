clear all;
% Point of this is to keep spacing of events but in non high-coflux points

bins = [1, 0.9, 0.75, 0.5, 0.25, 0.05];
subjects = {'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07','MSC09', 'MSC10'};
sessions = 1:10;
session_type = 'rest';
thresh = 0.05;
offset = -10:10;
compare = [];
total_compare = [];

% iterate over 100 times so we can average out the randomness of which 5
% points it will select 
for a = 1:100
    disp(a)
    for i = 1:length(bins)
        count = 0;
        for j = 1:length(subjects)
            for k = 1:length(sessions)      
                functionals_path = ['/projects/b1081/MSC/TaskFC/FC_Parcels/' session_type '/' subjects{j} '_parcel_timecourse.mat'];
                load(functionals_path)
                parcel_timeseries = parcel_time{sessions(k)};
                tmask = tmask_all{sessions(k)};
                parcel_timeseries(tmask==0, :) = NaN;
                timeseries = parcel_timeseries;
    
                numpts = round(thresh*size(timeseries(tmask==1, :),1));
                startpoint = floor((1-bins(i))*size(timeseries(tmask==1, :),1))+1;
                endpoint = startpoint + numpts -1;
                [keepfc, events_idx, timeseries_z] = getEventsFc(timeseries, startpoint, endpoint);
                rows = [];
    
                % for each subject/session/bin we are going to find events and
                % then offset them by 1-10 each direction and count how many
                % times this results in hitting a null point (stored in rows) 
                for l = 1:length(offset)
                    shifted_events = events_idx + offset(l);
                    over = shifted_events > size(timeseries, 1);
                    shifted_events(over) = shifted_events(over) - size(timeseries, 1);   
                    under = shifted_events < 1;
                    shifted_events(under) = shifted_events(under) + size(timeseries, 1);
                    shifted_events_ts = timeseries_z(shifted_events, :);
                    rows(l,:) = any(isnan(shifted_events_ts),2);
                end
    
                % we only use time points where none of the 21 offsets hit a
                % scrubbed value 
                clean_events = find(sum(rows) == 0);
    
    
                % semi arbitrary filter where we are only going to consider
                % subject/session/bin combos where we have at least 5 points to
                % use (could change this to be more... currently can use 53/90
                % sessions)
    
                % then we take a random 5 of those points in each session to
                % use as our "events." Yes there is randomness here that is not
                % really getting accounted for... 
                if length(clean_events) >=5
                    fc = corr(timeseries, 'rows', 'complete');
                    sample_events = events_idx(datasample(clean_events, 5, 'Replace', false));
                    for m=1:length(offset)
                        shifted_events = sample_events + offset(m);
                        over = shifted_events > size(timeseries, 1);
                        shifted_events(over) = shifted_events(over) - size(timeseries, 1);   
                        under = shifted_events < 1;
                        shifted_events(under) = shifted_events(under) + size(timeseries, 1);
                        shifted_events_ts = timeseries_z(shifted_events, :);
                        sample_fc = corr(shifted_events_ts);
                        compare(i, j, k, m) = corr(fc(:), sample_fc(:));
                    end
                else 
                    compare(i, j, k, :) = NaN; % else we set it to null 
                end 
            end
        end
    end
    compare_clean = compare;
    clean_sessions= squeeze(any(isnan(compare(:, :, :, 1)), 1));

    % for each iteration, we limit it ONLY to use the same
    % subjects/sessions for all bins
    for i = 1:size(clean_sessions,1)
        for j = 1:size(clean_sessions,2)
            if clean_sessions(i,j) == 1
                compare_clean(:, i, j, :) = NaN;
            end
        end
    end
    total_compare(a,:,:,:,:) = compare_clean;
end

save('offset_analysis.mat', 'total_compare')

%% function 
function [keepfc, events_idx, timeseries_z] = getEventsFc(timeseries, startpoint, endpoint)
    timeseries_z = (timeseries - nanmean(timeseries)) ./ nanstd(timeseries);
    
  
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

    eventCofluxRms = sqrt(nansum(coflux.^2));
    [~, rms_idx] = sort(eventCofluxRms, 'descend');
    events_idx = rms_idx(startpoint:endpoint);
    keepfc = corr(timeseries_z(rms_idx(startpoint:endpoint),:));
end
