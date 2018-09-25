% This is an m-file to test undecimated wavelets

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

%% Perform the Undecimated Wavelet Transform with Biorthogonal Wavelet

% Biorthogonal 3.1 Wavelet
%h = [-0.3535533906 1.0606601718 1.0606601718 -0.3535533906];
%g = [-0.1767766953 0.5303300859 -0.5303300859 0.1767766953];

% Biorthogonal 1.3 Wavelet
%h = [-0.08838834764831845 0.08838834764831845 0.7071067811865476 0.7071067811865476 0.08838834764831845 -0.08838834764831845];
%g = [0.0 0.0 0.7071067811865476 -0.7071067811865476 0.0 0.0];

% Biorthogonal 1.5 Wavelet
%h = [0.01657281518405971 -0.01657281518405971 -0.12153397801643787 0.12153397801643787 0.7071067811865476 0.7071067811865476 ...
%0.12153397801643787 -0.12153397801643787 -0.01657281518405971 0.01657281518405971];
%g = [0.0 0.0 0.0 0.0 -0.7071067811865476 0.7071067811865476 0.0 0.0 0.0 0.0]; 

% Haar Wavelet
%h = [sqrt(2)/2 sqrt(2)/2];
%g = [-sqrt(2)/2 sqrt(2)/2];

% Biorthogonal 3.1 Reconstruction Wavelet
h = [0.1767766952966369 0.5303300858899107 0.5303300858899107 0.1767766952966369];
g = [-0.3535533905932738 -1.0606601717798214 1.0606601717798214 0.3535533905932738];

% Perform one iteration of the transform
[MLow3,MHigh3] = UWT(Spectra(:,2),h,g);

% Perform Peak Detection
posneg = sign(MHigh3);
index = []; count = 1;
for i = 1:A-1
    if posneg(i) == -1 && posneg(i+1) == 1
        index(count,1) = i+1;
        count = count + 1;
    end
end

% Get the peak values from the original spectra
peaks = [];
for i = 1: length(index)
    peaks(i,:) = Spectra(index(i),:);
end

% Plot the spectra and the found peaks
figure;
plot(Spectra(:,1),Spectra(:,2));
hold on
scatter(peaks(:,1),peaks(:,2));
title('2,4,6-trimethylphenol and peak detection 1 iterations');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

%% Perfrom multiple iterations of the transform
[MLow4,MHigh4] = UWT1D(Spectra(:,2),h,g,400);

signdetect = [];
for i = 1:400
    posneg = sign(MHigh4(:,i));
    signdetect = [signdetect, posneg];
end

% This finds the peaks at each iteration of the transform
index = [];
for i = 400:-1:1
    count = 1;
    for j = 1:A-1
        if signdetect(j,i) == -1 && signdetect(j+1,i) == 1
            index(count,i) = j+1;
            count = count + 1;
        end
    end
    b(i) = count-1;
end

% Find the closest index values
Nindex = [];
count = 1;
for i = 1: b(400)
    col = 400;
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

%% Plot the spectra and the found peaks
figure;
plot(Spectra(:,1),Spectra(:,2));
hold on
scatter(Npeaks(:,1),Npeaks(:,2));
title('2,4,6-trimethylphenol and peak detection 50 iterations');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest
