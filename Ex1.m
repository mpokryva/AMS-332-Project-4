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
title("Ex 1Part 1 - Poisson spike trains")
%%%%%%%%% Part 3 %%%%%%%%%%

% Calculate average firing rate

f_count = firing_count(ISIs);
avg_f_rate = f_count / (N * T) %%%%% Average firing rate %%%%%%%%%%%%%%

%%%%%%%%% Part 4 %%%%%%%%%%

dt = 0.2;
bin_fr = zeros(1, T/dt); % Amount of firings in each bin.

for i = 1 : T / dt
   temp = ISIsCum((ISIsCum >= (i - 1) * dt) & (ISIsCum <= i * dt)); % Find number of firings in bin.
   bin_fr(i) = length(temp) / (N * dt); % Calculate bin firing rate.
end

mean_bin_fr = mean(bin_fr)


%%%%%%%%%%%%%%%%%%%%%%%% TODO: Plot part 4. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Part 5 %%%%%%%%%%
figure(3)
cv = coeff_var(ISIs);
plot(1:N, cv)
xlabel("Trial number")
ylabel("Coefficient of variability")
title("Ex 1 Part 5 - CV across 50 trials")
mean_cv = mean(cv) % Should be close to 1.

%%%%%%%%% Part 6 %%%%%%%%%%

ff = fano_factor(ISIs) % Should be close to 1.

%%%%%%%%% Part 7 %%%%%%%%%%

T = [0.5, 1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
for i = 1 : length(T)
    
end


function ff = fano_factor(ISIs)
    t_mean = firing_count(ISIs) / size(ISIs, 1);
    t_var = firing_var(ISIs);
    ff = t_var / t_mean;
end

function cv = coeff_var(ISIs)
    N = size(ISIs, 1);
    cv = zeros(1, N);
    for i = 1 : N
       cv(i) = nanstd(ISIs(i,:)) / nanmean(ISIs(i,:)); 
    end
end

function y = firing_count(ISIs)
    y = sum(sum(~isnan(ISIs)));
end

function y = firing_var(ISIs)
    counts = zeros(1, size(ISIs, 1));
    for i = 1 : size(ISIs, 1)
       c = firing_count(ISIs(i, :));
       counts(i) = c;
    end
    y = var(counts);
end


function [ISI, cumISI] = calculateISIs(lambda, T, max_size)
    ISI = zeros(1, max_size);
    for i = 1 : max_size
        ISI(i) = -log(rand(1,1))/lambda;
    end
    cumISI = cumsum(ISI);
    cumISI(cumISI > T) = NaN;
    ind = isnan(cumISI) == 1;
    ISI(ind) = NaN;
end