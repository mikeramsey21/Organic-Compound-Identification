function [Low,High] = UHWT(x)
% This is a matlab function programming the Undecimated Wavelet Transform

% Compute the length of the signal
Len = length(x);

% Compute the Undecimated HWT
Low(Len) = x(1) + x(Len);
High(Len) = x(1) - x(Len);
for i = 1:Len-1
    Low(i) = x(i) + x(i+1);
    High(i) = x(i) - x(i+1);
end

% Make the wavelet transfrom othgonal
Low = .5*Low';
High = .5*High';
end

