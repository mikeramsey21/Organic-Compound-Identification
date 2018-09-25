function [Low, High] = UWT(x,h,g)
% This is a matlab function that performs the othogonal wavelet transform
% with an orthogonal filter h of even length of a vector x

% This function inputs a vecor x and an orhogonal filter h
% This function outpus the Lowpass and Highpass portions of the Undecimated
% Wavelet transform for vector x

% Note that this function only performs one iteration of the transform

% If we have a row vector, tranform to a column vector
[A,B] = size(x);
if A == 1;
    x = x';
end

% Compute the length of the signal and the filter
Vectorlen = length(x);
Filterlen = length(h);

% Compute the number of wrapping rows
wrapnum = Filterlen - 1;

% Compute the wavelet transform for the non-wrapping rows
for i = 1: Vectorlen - wrapnum
    Low(i) = dot(x(i:i+wrapnum),h);
    High(i) = dot(x(i:i+wrapnum),g);
end

% Compute the wavelet tranform for the wrapping rows
for i = Vectorlen-wrapnum + 1: Vectorlen
    temp = [x(i:Vectorlen); x(1: Filterlen - (Vectorlen - i + 1))];
    Low(i) = dot(temp,h);
    High(i) = dot(temp,g);
end

% Transpose the vectors into column vectors
Low = Low';
High = High';
    

end

