function [startQ goalQ] = generateRandomConf()
    startQ = rand(1,5)*2*pi;
    goalQ = rand(1,5)*2*pi;
end