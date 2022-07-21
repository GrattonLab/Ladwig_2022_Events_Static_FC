clear all;
load('PATH TO SIM DATA YOU CREATED')


sessions = 1:10;
subjects = {'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07','MSC09', 'MSC10'};
thresh = 0.05;
bins = [1, thresh];

all_events = [];
non_events = [];

for i = 1:length(bins)
    for k = 1:length(subjects)
        subject = subjects{k};
        for j = 1:length(sessions)
            bin = bins(i);
            timeseries =fake_signals{k,j}';
            if length(timeseries) == 0
                break;
            end
            numpts = round(thresh*size(timeseries,1));
            startpoint = floor((1-bin)*size(timeseries,1))+1;
            endpoint = startpoint + numpts-1;
            [keepfc, events_idx, rms_idx, timeseries_z, eventCofluxRms] = getEventsFc(timeseries, startpoint, endpoint);
            events = timeseries_z(events_idx,:);
            if(i == 1)
                    all_events = [all_events; events];
            else 
                    non_events = [non_events; events];
            end
        end
    end
end

%% get PC1
all_pts = [all_events non_events];
niche = corr(all_pts);
[coeff,score,latent, tsquared, explained] = pca(niche);
pc1 = coeff(:,1);
pc1_sim = pc1;

figure()
scatter(ones(1,333), abs(pc1(1:333)), 'jitter','on','jitterAmount',0.15, 'MarkerFaceAlpha',0.3')
hold on 
boxplot(abs(pc1(1:333)))
scatter(2*ones(1,333), abs(pc1(334:666)), 'jitter','on','jitterAmount',0.15, 'MarkerFaceAlpha',0.3')
boxplot(abs(pc1(334:666)), 'positions', 2)
xlim([0 3])

save('pca_sim.mat', 'pc1_sim')


%% load a cifti, match up 333 with vertics, plot it 

cifti = ft_read_cifti_mod('PATH TO A CIFTI');
cifti_copy = cifti;
for i = 1:333
    cifti_copy.data(find(cifti.data==i)) = pc1(i);
end

ft_write_cifti_mod('pc1_sim.dtseries.nii',cifti_copy);

%%  events code 
function [keepfc, events_idx, rms_idx, timeseries_z, eventCofluxRms] = getEventsFc(timeseries, startpoint, endpoint)
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