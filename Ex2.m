close all

f_max = 200; % 200 spikes/sec.
v_spk = 20; % mV
N = 100; % Number of neurons in network.
tau = 10 * 10 .^ -3; % ms
max_f = 150;
df = 0.01;
equilJ = zeros(2, max_f/df); 
phi_f = zeros(1, max_f/df);
f = (0 : df : max_f);
J = [0.1, 0.185, 0.3, 0.4];
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

%%%%%%%%% Part 2 %%%%%%%%%%%%%
dt = 0.1 * 10 .^ -3; % ms
tau_f = 10 * 10 .^ -3 % ms
tau = 10 * 10 .^ -3; % ms
t = 500 * 10 .^ -3 % ms
f_0 = [5, 55, 75, 130];
J = [0.1, 0.185, 0.3, 0.4];
%f_0 = [25];
time = zeros(1, t/dt);
f = zeros(1, t/dt);
for k = 1 : length(J)
    for i = 1 : length(f_0)
        figure(i+1);
        hold on;
        f(1) = f_0(i);
        time(1) = 0;
        for j = 1 : t/dt
            mv = mu_v(J(k), N, f(j), tau);
            sv = sigma_v(J(k), N, f(j), tau); 
            resp = response(mv, sv, v_spk, f_max, N);
            df = (((f(j) * -1) + resp) / tau_f) * dt;
            f(j+1) = f(j) + df;
            time(j+1) = time(j) + dt;
        end
        plot(time, f);
        xlabel("Time (s)")
        ylabel("Firing rate (spk/s)")
        title("Firing rate with J = " + J(k) + ", f_0 = " + f_0(i));
    end 
end

%legend_labels = strings(1, length(f_0));
%for i = 1 : length(legend_labels)
 %  legend_labels(i) = strcat("f_0=", num2str(f_0(i)));
%end
%lgd = legend(legend_labels);
%lgd.FontWeight = "bold";



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