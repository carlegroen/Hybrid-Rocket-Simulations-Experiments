clc
close all
clear all
format long g

T_amb = 273 + 20;
T_end = 273 + 2225;
R = 8.3145;
M = 1;
gamma = 1.2;
T = linspace(T_amb,T_end,1000)
v = sqrt(R*T/M * gamma)
plot(T,v)

gas_density = 0.1

n_dot_in = 1.71609682164548;

m_dot =  gas_density * v *  0.002*0.002;


figure(2)
plot(T,m_dot)

n_dot_out = m_dot/5;
n_dot_in  = 1.1*n_dot_out;
n = n_dot_in-n_dot_out;

At = 0.002*0.002
P = m_dot .* sqrt(T) .* 1/(At * sqrt(gamma/R) .* ((gamma+1)/2).^(-(gamma+1)/(2.*(gamma-1))))

figure(3)
plot(T,P)

plot((T/T(end)),((T/T(end)).^(-1)).^(-gamma/(gamma-1)))