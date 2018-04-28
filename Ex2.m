close all

f_max = 200; % 200 spikes/sec.
v_spk = 20; % mV
N = 100; % Number of neurons in network.
tau = 10 .^ -3; % ms

max_size = 100;
equilJ = zeros(1, max_size);
i = 1;
for f = 1 : 150
    for J = 2 : 1000
       mv = mu_v(J, N, f, tau);
       sv = sigma_v(J, N, f, tau);
       res = response(mv, sv, v_spk, f_max, N);
       if res == f
           equilJ(i) = J;
           i = i + 1;
       end
    end  
end


function y = mu_v(J, N, f, tau)
    y = J * N * f * tau;
end

function y = sigma_v(J, N, f, tau)
    y = J * sqrt(N * f * tau);
end

function y = response(mu_v, sigma_v, v_spk, f_max, N)
    num = f_max;
    ex = (-sqrt(2) * (mu_v - v_spk))/(sigma_v * sqrt(N));
    den = 1 + exp(ex);
    y = num / den;
end