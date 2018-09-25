function [ Npeaks ] = peaks(Spectra,N)
% This is a matlab function that performs peak detection using the
% Undecimated HWT

% Perform the Undecimated Transform
[MLow2,MHigh2] = UHWT1D(Spectra(:,2),N);
[A,B] = size(Spectra);

signdetect = [];
for i = 1:N
    posneg = sign(MHigh2(:,i));
    signdetect = [signdetect, posneg];
end

% This finds the peaks at each iteration of the transform
index = [];
for i = N:-1:1
    count = 1;
    for j = 1:A-1
        if signdetect(j,i) == 1 && signdetect(j+1,i) == -1
            index(count,i) = j+1;
            count = count + 1;
        end
    end
    b(i) = count-1;
end

% Find the closest index values
Nindex = [];
count = 1;
for i = 1: b(N)
    col = N;
    temp = index(i,col);
    while col ~= 1
        sub = abs(temp - index(:,col-1));
        minimum = min(sub);
        new = find(sub==minimum);
        temp = index(new,col-1);
        col = col - 1;
    end
    Nindex(count) = temp;
    count = count + 1;
end

% Get the peak values from the original spectra
for i = 1: length(Nindex)
    Npeaks(i,:) = Spectra(Nindex(i),:);
end

% Delete the last row which is somehow redundant
[A,B] = size(Npeaks);
Npeaks(A,:) = [];
end

