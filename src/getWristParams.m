close all
clear, clc

U = zeros(1,16);
H = zeros(1,16);
for i = 1:4
    for j = 1:4
        [w, H((i-1)*4+j), U((i-1)*4+j)] = notchsynthesis(0.005 * j, 0.003 * i);
    end
end
    