close all
clear
clc

%% Test Skript

x = zeros(10, 1);
y = zeros(10, 1);
A = zeros(10);
alpha = 3;
beta  = 5;

for i = 0:9
    x(i+1) = i;
    y(i+1) = i*2;
    for j = 0:9
        A(i+1, j+1) = i * 10 + j;
    end
end

y = beta * y + alpha * A * x 