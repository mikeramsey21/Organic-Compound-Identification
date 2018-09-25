function [Low,High] = UHWT1D(x,its)
% This is a matlab function programming multiple iterations of the
% Undecimated Haar Wavelet Transform

High = []; Low = [];
for j = 1:its
    [Low1,High1] = UHWT(x);
    x = Low1;
    High = [High, High1];
    Low = [Low, Low1];
end

end

