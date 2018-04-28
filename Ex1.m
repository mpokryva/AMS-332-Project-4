close all

%%%%%%%%% Part 1 %%%%%%%%%%
N = 50;
lambda = 10; % Spikes/sec
T = 5;
max_size = 100;
ISIsCum = zeros(N, max_size);
ISIs = zeros(N, max_size);
for i = 1 : N
    [ISIs(i,:), ISIsCum(i,:)] = calculateISIs(lambda, T, max_size);
end

%%%%%%%%% Part 2 %%%%%%%%%%

figure(1)
plot(ISIsCum, 1:N, '.k')
xlabel("Time (s)")
ylabel("Trial number")
title("Part 1 - Poisson spike trains")
%%%%%%%%% Part 3 %%%%%%%%%%

% Calculate average firing rate

nans = isnan(ISIsCum);
nan_occ = sum(sum(nans)); % Number of NaNs
total_spikes = (N * max_size) - nan_occ;
avg_f_rate = total_spikes / (N * T) %%%%% Average firing rate %%%%%%%%%%%%%%

%%%%%%%%% Part 4 %%%%%%%%%%

dt = 0.2;
bin_fr = zeros(1, T/dt); % Amount of firings in each bin.

for i = 1 : T / dt
   temp = ISIsCum((ISIsCum >= (i - 1) * dt) & (ISIsCum <= i * dt)); % Find number of firings in bin.
   bin_fr(i) = length(temp) / (N * dt); % Calculate bin firing rate.
end

mean_bin_fr = mean(bin_fr)

%% TODO: Plot part 4.

%%%%%%%%% Part 5 %%%%%%%%%%

cv = zeros(1, N);
for i = 1 : N
   cv(i) = nanstd(ISIs(i,:)) / nanmean(ISIs(i,:)); 
end

cv_mean = mean(cv) % Should be close to 1.



function [ISI, cumISI] = calculateISIs(lambda, T, max_size)
    ISI = zeros(1, max_size);
    for i = 1 : max_size
        ISI(i) = -log(rand(1,1))/lambda;
    end
    cumISI = cumsum(ISI);
    cumISI(cumISI > T) = NaN;
    ind = find(isnan(cumISI) == 1);
    ISI(ind) = NaN;
end