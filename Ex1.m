close all

N = 50;
lambda = 10; % Spikes/sec
T = 5;
max_size = 100;
ISIs = zeros(N, max_size);
for i = 1 : N
    ISIs(i,:) = calculateISI(lambda, T, max_size);
end
disp(ISIs)



function ISI = calculateISI(lambda, T, max_size)
    ISI = zeros(1, max_size);
    for i = 1 : max_size
        ISI(i) = -log(rand(1,1))/lambda;
    end
    ISI = cumsum(ISI);
    ISI(ISI > T) = NaN;
end