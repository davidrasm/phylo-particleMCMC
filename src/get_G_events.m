function [G_lineages, G_coal_events, G_lineages_dt, G_indices, G_dt_ref, event_times] = get_G_events(MCMC_params, sample_times, sample_sizes, coal_times, obsv_times)
%set up data structures for genealogical data
%major update 100410: data is now structured by dt interval not events

%All times should be in months since tEnd
sample_times = sort(sample_times);
coal_times = coal_times(coal_times < max(obsv_times)); %remove coal times that occur before t0
coal_times = sort(coal_times);
obsv_times = sort(obsv_times);
lineages = 0;
tEnd = max(obsv_times);
dt_times = 0:MCMC_params.dt:tEnd;

event_times = [coal_times, obsv_times, sample_times, dt_times];
event_times = sort(unique(event_times));

G_lineages = zeros(1, length(event_times));
G_coal_events = zeros(1, length(event_times));
G_dt_ref = zeros(1, length(event_times));
n = 0;
for i = 1:length(event_times)
    
    time = event_times(i);
    
    %Update G_lineages
    if ismember(time, sample_times)
        loc = find(sample_times == time);
        lineages = lineages + sample_sizes(loc);
        %G_lineages(i) = lineages;
    end

    %Update G_coal_events
    if ismember(time, coal_times)
        G_coal_events(i) = 1;
        coal_locs = find(coal_times == time); %modified to handle more than one coal event at one time
        lineages = lineages - length(coal_locs);
    end
    
    if ismember(time, obsv_times)
        n = n + 1;
        G_indices(n) = i;
    end

    G_lineages(i) = lineages;
    
    %%%Find closest dt time%%%%
    diff = time - dt_times;
    abs_diff = abs(diff);
    dt_loc = find(abs_diff == min(abs_diff));
    G_dt_ref(i) = dt_loc; 
    
end


G_lineages_dt = zeros(1, length(dt_times));
for j = 1:length(dt_times) 
    time = dt_times(j);
    loc = find(event_times == time);
    n = n + 1;
    G_lineages_dt(j) = G_lineages(loc);
end
    
%%%Convert times to times since t0
G_lineages = fliplr(G_lineages);
G_coal_events = fliplr(G_coal_events);
G_lineages_dt = fliplr(G_lineages_dt);
G_indices = sort(length(event_times) - G_indices + 1);
G_dt_ref = sort(length(dt_times) - G_dt_ref + 1);
event_times = sort(tEnd - event_times);


