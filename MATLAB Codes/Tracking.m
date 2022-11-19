%% Author: V. Ruban Vishnu Pandian, EE19B138
%% EE6110: Adaptive Signal Processing - Project
%% Date: 17/11/2022

%% Standard LMS

a = sqrt(1.3)/sqrt(3);
b = [a 1/sqrt(2) 1/sqrt(6); a -1/sqrt(2) 1/sqrt(6); a 0 -2/sqrt(6)];
b = b/100;

N = 2500;
a = [1; 0.6; -0.1; 0.1; 0.1];
sigv = 0.5;

V_vec = sqrt(sigv)*randn(1,N+2);
num = 1;
denom = a;

Z_vec = filter(num,denom,V_vec);

u = 3.402*1e-2;

sigU = 0.1;
U1 = exprnd(sqrt(0.5*sigU),1,N);
U2 = exprnd(sqrt(0.5*sigU),1,N);
U = U1-U2;

A = b; 
p = randn(3,N);
q = A*p;

e = zeros(1,N);
h = zeros(3,N+1);
w = zeros(3,N+1);
MSD = zeros(1,N);

h(:,1) = [-1; 0; 1];

for i=1:N
    curr_Y = (fliplr(Z_vec(i:i+2))).';
    curr_x = (h(:,i)'*curr_Y)+U(i);

    e(i) = curr_x-(w(:,i)'*curr_Y);
    w(:,i+1) = w(:,i)+(u*e(i)'*curr_Y);
    h(:,i+1) = h(:,i)+q(:,i);
    MSD(i) = norm(w(:,i+1)-h(:,i+1))^2;
end

figure(1)
plot(1:N,h(1,2:end),'blue');
hold on;
plot(1:N,h(2,2:end),'blue');
plot(1:N,h(3,2:end),'blue');
scatter(1:N,w(1,2:end),1,'filled','red');
scatter(1:N,w(2,2:end),1,'filled','red');
scatter(1:N,w(3,2:end),1,'filled','red');
legend('True system','','','LMS tracking');
xlabel('Time index');
ylabel('Weights amplitude');
title("Weights variation as function of time for LMS");

figure(2)
plot(1:N,10*log10(MSD));
xlabel('Time index');
ylabel('MSD');
title("MSD variation as function of time for LMS");

%% With error weighing function
a = sqrt(1.3)/sqrt(3);
b = [a 1/sqrt(2) 1/sqrt(6); a -1/sqrt(2) 1/sqrt(6); a 0 -2/sqrt(6)];
b = b/100;

N = 2500;
a = [1; 0.6; -0.1; 0.1; 0.1];

u = 3.402*1e-2;

sigv = 0.5;
V_vec = sqrt(sigv)*randn(1,N+2);

num = 1;
denom = a;

Z_vec = filter(num,denom,V_vec);

sigU = 0.1;
U1 = exprnd(sqrt(0.5*sigU),1,N);
U2 = exprnd(sqrt(0.5*sigU),1,N);
U = U1-U2;

A = b; 
p = randn(3,N);
q = A*p;

e = zeros(1,N);
h = zeros(3,N+1);
w = zeros(3,N+1);
MSD = zeros(1,N);

h(:,1) = [-1; 0; 1];

p = 3;
ep = 1e-4;

for i=1:N
    curr_Y = (fliplr(Z_vec(i:i+2))).';
    curr_x = (h(:,i)'*curr_Y)+U(i);

    e(i) = curr_x-(w(:,i)'*curr_Y);

    f = (abs(e(i))^p)/(ep+(abs(e(i)))^(p+1));
    w(:,i+1) = w(:,i)+(u*f*e(i)'*curr_Y);
    h(:,i+1) = h(:,i)+q(:,i);
    MSD(i) = norm(w(:,i+1)-h(:,i+1))^2;
end

figure(3)
plot(1:N,h(1,2:end),'blue');
hold on;
plot(1:N,h(2,2:end),'blue');
plot(1:N,h(3,2:end),'blue');
scatter(1:N,w(1,2:end),1,'filled','red');
scatter(1:N,w(2,2:end),1,'filled','red');
scatter(1:N,w(3,2:end),1,'filled','red');
legend('True system','','','LMS tracking');
xlabel('Time index');
ylabel('Weights amplitude');
title("Weights variation as function of time for NFRMS");

figure(4)
plot(1:N,10*log10(MSD));
xlabel('Time index');
ylabel('MSD');
title("MSD variation as function of time for NFRMS");