%% Code to generate static simulation for all MSC subjects as in Fig 2
clear all; 
subjects = {'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07','MSC09', 'MSC10'};
sessions = 1:10;
session_type = 'rest';
thresh = 0.05;
bins = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, thresh];
mod_threshold = 0.05;
fake_signals = {};

%% generate fake data for each subject, calculate mod and fc 
for i = 1:length(bins)
    disp(i) 
    count = 0;
    for j = 1:length(subjects)
        for k = 1:length(sessions)

            bin = bins(i);
            subject = subjects{j};
            session = sessions(k);

            functionals_path = ['YOUR FILE PATH' session_type '/' subject '_parcel_timecourse.mat'];
            load(functionals_path)
            parcel_timeseries = parcel_time{session};
            tmask = tmask_all{session};
            timeseries = parcel_timeseries(tmask==1, :)';

            if(size(timeseries,2) <= 333)
                fc_sim(i,j,k) = NaN;
                mod_sim(i,j,k) = NaN;
                continue;
            end

            real_signal = timeseries;
            real_fc = corr(real_signal');
            rand_signal = randn(size(real_signal));
            
            
            % make sure that the matrix is positive-determinant
            try
               R = chol(real_fc);
            catch
               fc_sim(i,j, k) = NaN;
               mod_sim(i,j,k) = NaN;
               continue;
            end  
            
            % get eigenvectors of real FC 
            [V D] = eig(real_fc);
            V_use = V*sqrt(D);

            % apply them to our random signal calculate coflux, modularity
            % etc
            fake_signal = V_use*rand_signal;
            fake_signals{j,k} = fake_signal;
            fake_fc = corr(fake_signal');

            numpts = round(thresh*size(fake_signal,2));
            startpoint = floor((1-bins(i))*size(fake_signal,2))+1;
            endpoint = startpoint + numpts -1;

            keepfc = getEventsFc(fake_signal', startpoint, endpoint);
            fc_sim(i, j, k) = corr(fake_fc(:), keepfc(:));
            
            [mat_thresh r kden] = matrix_thresholder(keepfc, mod_threshold ,'kden');
            mat_thresh(find(mat_thresh > 0))=1;
            [Ci Q] = modularity_und(mat_thresh);
            mod_sim(i, j, k) = Q;
        end
    end
end

% save out our sim data the same way we did the real 
save('sim_output.mat', 'mod_sim', 'fc_sim', 'fake_signals')

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

