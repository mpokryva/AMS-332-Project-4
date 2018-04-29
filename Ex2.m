close all

f_max = 200; % 200 spikes/sec.
v_spk = 20; % mV
N = 100; % Number of neurons in network.
tau = 10 .^ -3; % ms
max_f = 150;
df = 0.01;
equilJ = zeros(2, max_f/df); 
phi_f = zeros(1, max_f/df);
f = (0 : df : max_f);
J = [0.5, 0.815, 1, 1.5, 2];
%J = (0 : 0.1 : 3);
k = 1;
hold on
for i = 1 : length(J)
    for j = 1 : length(f)
       mv = mu_v(J(i), N, f(j), tau);
       sv = sigma_v(J(i), N, f(j), tau);    
       phi_f(j) = response(mv, sv, v_spk, f_max, N);
    end
    plot(f, phi_f);
end
plot(f, f, "--");
xlabel("f (spk/s)")
ylabel("\phi(f) (spk/s)")
title("Ex 2 Part 1 - \phi(f) vs. f")

legend_labels = strings(1, length(J));
for i = 1 : length(legend_labels)
   legend_labels(i) = strcat("J=", num2str(J(i))); 
end
lgd = legend(legend_labels);
lgd.FontWeight = "bold";
hold off


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