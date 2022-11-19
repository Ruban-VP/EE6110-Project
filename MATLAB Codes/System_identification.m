%% Author: V. Ruban Vishnu Pandian, EE19B138
%% EE6110: Adaptive Signal Processing - Project
%% Date: 18/11/2022

clear;
close all;

%% System identification
N = 1e4;
f = [1; 2; 3];
L = length(f);

sigy = 0.5;
y = sqrt(sigy)*randn(1,N+L-1);
x = conv(f,y);
x = x(L:N+L-1);

sigU = 5;
U1 = exprnd(sqrt(0.5*sigU),1,N);
U2 = exprnd(sqrt(0.5*sigU),1,N);
U = U1-U2;

x = x+U;

w = zeros(length(f),1);
u = 0.1;

pow = 2;
ep = 1e-4;

MSD = zeros(1,N);

for i=1:N
    Y = (fliplr(y(i:i+L-1))).';
    X = x(i);
    e = X-(w'*Y);

    fun = (abs(e)^pow)/(ep+abs(e)^(pow+1));
    %fun = 1;
    w = w + u*fun*e'*Y;
    MSD(i) = norm(f-w)^2;
end

plot(1:N,MSD);
