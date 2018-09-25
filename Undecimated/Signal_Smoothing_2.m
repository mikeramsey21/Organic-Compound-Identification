% This is an m-file to perfrom signal smoothing

% Reading in our IR Spectra
% The first column is the associated wavelength for each transmittance. The
% second column is the % transmittance 
Spectra = csvread('2-4-6-trimethylphenol.CSV');

% There is alot of useless information in these files. Therefore we are
% going to trim the beginning of the data set
[A,B] = size(Spectra);
while Spectra(1,2) == 0
    Spectra(1,:) = [];
end
[A,B] = size(Spectra);
Spectra(A,:) = [];

% Re-compute the size of Spectra and generate a random vecotor for noise
% incorporation
[A,B] = size(Spectra);

% Plot the original Spectra for 2-4-6-trimethylphenol and Detail
% Coefficients
figure;
plot(Spectra(:,1),Spectra(:,2));
hold on
title('2,4,6-trimethylphenol');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

%% Perform signal smoothing
% The divisors of 7208
% 1, 2, 4, 8, 17, 34, 53, 68, 106, 136, 212, 424, 901, 1802, 3604, 7208
piece = 424;
smooth = ones(1,piece)/piece;
smoothsignal = filter(smooth,1,Spectra(:,2));

% Plot the original and smoothed signal
figure;
plot(Spectra(:,1),[smoothsignal Spectra(:,2)]);
title('2,4,6-trimethylphenol');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

% Adjust for filter delay
fdelay = (piece-1)/4;
figure;
plot(Spectra(:,1),Spectra(:,2));
hold on
plot(Spectra(:,1)-fdelay, smoothsignal);
title('2,4,6-trimethylphenol');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

% Update the smoothed Spectra
Spectra2 = [Spectra(:,1)-fdelay, smoothsignal];

% Perform a a second iteration of signal smoothing
piece2 = 106;
smooth = ones(1,piece2)/piece2;
smoothsignal = filter(smooth,1,Spectra2(:,2));

% Adjust for filter delay
fdelay2 = (piece2-1)/4;
figure;
plot(Spectra(:,1),Spectra(:,2));
hold on
plot(Spectra2(:,1)-fdelay2, smoothsignal);
title('2,4,6-trimethylphenol');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

% Update the smoothed Spectra
Spectra3 = [Spectra2(:,1)-fdelay2, smoothsignal];

%% Perform peak detection using UHWT
% Perform n iterations of the Undecimated Transform and plot the result
[MLow2,MHigh2] = UHWT1D(Spectra3(:,2),100);

signdetect = [];
for i = 1:100
    posneg = sign(MHigh2(:,i));
    signdetect = [signdetect, posneg];
end

% This finds the peaks at each iteration of the transform
index = [];
for i = 100:-1:1
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
for i = 1: b(100)
    col = 100;
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
    Npeaks(i,:) = Spectra3(Nindex(i),:);
end        

% Plot the spectra and the found peaks
figure;
plot(Spectra3(:,1),Spectra3(:,2));
hold on
scatter(Npeaks(:,1),Npeaks(:,2));
hold on
plot(Spectra(:,1),Spectra(:,2));
title('2,4,6-trimethylphenol and peak detection 100 iterations');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest
