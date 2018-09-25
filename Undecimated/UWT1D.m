function [Low,High] = UWT1D(x,h,g,its)
% This is a matlab function programming multiple iterations of the
% Undecimated Wavelet Transform with associated filter h

High = []; Low = [];
for j = 1:its
    [Low1,High1] = UWT(x,h,g);
    x = Low1;
    High = [High, High1];
    Low = [Low, Low1];
end

end