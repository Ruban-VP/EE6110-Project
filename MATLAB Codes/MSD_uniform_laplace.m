%% Author: V. Ruban Vishnu Pandian, EE19B138
%% EE6110: Adaptive Signal Processing - Project
%% Date: 18/11/2022

clear;
close all;

%% System identification for 20% uniform and 80% laplace

N = 1e3;
N_trials = 1e3;
L = 5;
f = randn(1,L);
L = length(f);

sigy = 1;
y = sqrt(sigy)*randn(1,N+L-1);
x = conv(f,y);
x = x(L:N+L-1);

MSD_1 = zeros(1,N);
MSD_2 = zeros(1,N);
MSD_3 = zeros(1,N);

pow = 2;
ep = 1;

u_1 = 0.0004;
u_2 = 0.0012;
u_3 = 0.001;

for ind=1:N_trials
    
    w_1 = zeros(length(f),1);
    w_2 = zeros(length(f),1);
    w_3 = zeros(length(f),1);

    sigU = 10;
    unif_noise = (rand(1,N)*2-1)*sqrt(3); 
    U1 = exprnd(sqrt(0.5*sigU),1,N);
    U2 = exprnd(sqrt(0.5*sigU),1,N);
    lap_noise = U1-U2;
    
    noise = 0.2*unif_noise+0.8*lap_noise; 
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

%% System identification for 50% uniform and 50% laplace

N = 1e3;
N_trials = 1e3;
L = 5;
f = randn(1,L);
L = length(f);

sigy = 1;
y = sqrt(sigy)*randn(1,N+L-1);
x = conv(f,y);
x = x(L:N+L-1);

MSD_1 = zeros(1,N);
MSD_2 = zeros(1,N);
MSD_3 = zeros(1,N);

pow = 2;
ep = 1;

u_1 = 0.0004;
u_2 = 0.0012;
u_3 = 0.001;

for ind=1:N_trials
    
    w_1 = zeros(length(f),1);
    w_2 = zeros(length(f),1);
    w_3 = zeros(length(f),1);

    sigU = 10;
    unif_noise = (rand(1,N)*2-1)*sqrt(3); 
    U1 = exprnd(sqrt(0.5*sigU),1,N);
    U2 = exprnd(sqrt(0.5*sigU),1,N);
    lap_noise = U1-U2;
    
    noise = 0.5*unif_noise+0.5*lap_noise; 
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

figure(2)
plot(1:N,10*log10(MSD_1));
hold on;
plot(1:N,10*log10(MSD_2));
plot(1:N,10*log10(MSD_3));
legend('LMLS','LLAD','NFRMS');

%% System identification for 80% uniform and 20% laplace

N = 1e3;
N_trials = 1e3;
L = 5;
f = randn(1,L);
L = length(f);

sigy = 1;
y = sqrt(sigy)*randn(1,N+L-1);
x = conv(f,y);
x = x(L:N+L-1);

MSD_1 = zeros(1,N);
MSD_2 = zeros(1,N);
MSD_3 = zeros(1,N);

pow = 2;
ep = 1;

u_1 = 0.0004;
u_2 = 0.0012;
u_3 = 0.001;

for ind=1:N_trials
    
    w_1 = zeros(length(f),1);
    w_2 = zeros(length(f),1);
    w_3 = zeros(length(f),1);

    sigU = 10;
    unif_noise = (rand(1,N)*2-1)*sqrt(3); 
    U1 = exprnd(sqrt(0.5*sigU),1,N);
    U2 = exprnd(sqrt(0.5*sigU),1,N);
    lap_noise = U1-U2;
    
    noise = 0.8*unif_noise+0.2*lap_noise; 
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

figure(3)
plot(1:N,10*log10(MSD_1));
hold on;
plot(1:N,10*log10(MSD_2));
plot(1:N,10*log10(MSD_3));
legend('LMLS','LLAD','NFRMS');