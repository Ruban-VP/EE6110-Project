%% Author: V. Ruban Vishnu Pandian, EE19B138
%% EE6110: Adaptive Signal Processing - Project
%% Date: 19/11/2022

clear;
close all;

%% System identification for 20% uniform and 80% laplace

N = 1e3;
N_trials = 1e2;
L = 16;
f = randn(L,1);
f = f/sqrt(f'*f);
L = length(f);

sigy = 1;
y = sqrt(sigy)*randn(1,N+L-1);

inp_fil_num = 1;
inp_fil_den = [1; -0.8];
y = filter(inp_fil_num,inp_fil_den,y);

x = conv(f,y);
x = x(L:N+L-1);

MSD_1 = zeros(1,N);
MSD_2 = zeros(1,N);
MSD_3 = zeros(1,N);

pow = 2;
ep = 0.1;

u_1 = 0.005;
u_2 = 0.003;
u_3 = 0.0047;

for ind=1:N_trials
    
    w_1 = zeros(length(f),1);
    w_2 = zeros(length(f),1);
    w_3 = zeros(length(f),1);

    bino = binornd(1,0.005,1,N);
    gauss = sqrt(5000)*randn(1,N);
    noise = bino.*gauss;

    x_pow = (x*x')/N;
    SNR = 10^(20/10);
    noise_pow = x_pow/SNR;

    x = x+(2*rand(1,N)-1)*sqrt(3*noise_pow);
    x = x+noise;

    for i=1:N
        Y = (fliplr(y(i:i+L-1))).';
        X = x(i);
        e1 = X-(w_1'*Y);
        e2 = X-(w_2'*Y);
        e3 = X-(w_3'*Y);
    
        % LMLS
        fun = (abs(e1)^2)/(ep+abs(e1)^2);
        w_1 = w_1 + u_1*fun*e1'*Y;
        MSD_1(i) = MSD_1(i)+norm(f-w_1)^2;
    
        % LLAD
        fun = 1/(ep+abs(e2));
        w_2 = w_2 + u_2*fun*e2'*Y;
        MSD_2(i) = MSD_2(i)+norm(f-w_2)^2;
    
        % NFRMS
        fun = (abs(e3)^pow)/(ep+abs(e3)^(pow+1));
        w_3 = w_3 + u_3*fun*e3'*Y;
        MSD_3(i) = MSD_3(i)+norm(f-w_3)^2;
    end
end

MSD_1 = MSD_1/N_trials;
MSD_2 = MSD_2/N_trials;
MSD_3 = MSD_3/N_trials;

figure(1)
plot(1:N,10*log10(MSD_1));
hold on;
plot(1:N,10*log10(MSD_2));
plot(1:N,10*log10(MSD_3));
legend('LMLS','LLAD','NFRMS');