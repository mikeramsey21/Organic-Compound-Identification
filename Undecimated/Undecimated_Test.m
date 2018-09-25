% This is an m-file to test undecimated wavelets

% Reading in our IR Spectra
% The first column is the associated wavelength for each transmittance. The
% second column is the % transmittance 
Spectra = csvread('4-ethylphenol.CSV');

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

[Low,High] = UHWT(Spectra(:,2));

% Plot the original Spectra for 2-4-6-trimethylphenol and Detail
% Coefficients
figure;
plot(Spectra(:,1),Spectra(:,2),'LineWidth',1.5);
hold on
plot(Spectra(:,1),High,'LineWidth',1.5);
title('4-ethylphenol');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest
legend('Original Spectra','Detail Coefficients of 1st Transform');
xlabel('Wavenumber (cm^-^1)');
ylabel('% Absorbtion');
% We can see that zero crossings of the detail coefficients generally
% correspond to peaks in the spectra

% Perform two iterations of the Undecimated Transform and plot the result
[MLow,MHigh] = UHWT1D(Spectra(:,2),4);
figure;
plot(Spectra(:,1),MLow(:,3),'LineWidth',1.5);
hold on
plot(Spectra(:,1),MHigh(:,4),'LineWidth',1.5);
title('4-ethylphenol');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest
legend('Averaging Coefficients of 3rd Transform','Detail Coefficients of 4th Transform');
xlabel('Wavenumber (cm^-^1)');
ylabel('% Absorbtion');

% Perform 100 iterations of the Undecimated Transform and plot the result
[MLow2,MHigh2] = UHWT1D(Spectra(:,2),300);
figure;
plot(Spectra(:,1),Spectra(:,2),'LineWidth',1.2);
hold on
plot(Spectra(:,1),MHigh2(:,10));
title('2,4,6-trimethylphenol and 100 iterations');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

%% We begin the naive peak detection algorithm
% A is the length of the spectra

% We first perform this peak detection on one iteration of the tranform
posneg = sign(High);
index = []; count = 1;
for i = 1:A-1
    if posneg(i) == 1 && posneg(i+1) == -1
        index(count) = i+1;
        count = count + 1;
    end
end

% Get the peak values from the original spectra
for i = 1: length(index)
    peaks(i,:) = Spectra(index(i),:);
end

% Plot the spectra and the found peaks
figure;
plot(Spectra(:,1),Spectra(:,2),'LineWidth',1.5);
hold on
scatter(peaks(:,1),peaks(:,2),'LineWidth',1.5);
title('4-ethylphenol with Peak Detection - 1 Iteration');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

%% Now perform peak detection algorithm on muliple iterations
signdetect = [];
for i = 1:10
    posneg = sign(MHigh(:,i));
    signdetect = [signdetect, posneg];
end

% This finds the peaks at each iteration of the transform
index = [];
for i = 10:-1:1
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
for i = 1: b(10)
    col = 10;
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

% Plot the spectra and the found peaks
figure;
plot(Spectra(:,1),Spectra(:,2),'LineWidth',1.5);
hold on
scatter(Npeaks(:,1),Npeaks(:,2),'LineWidth',1.5);
title('4-ethylphenol with Peak Detection - 10 Iterations');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

%% We proceed to using the daubechies wavelets for peak detection
%% Daubechies Wavelets Perform Poorly Unfortunately

% Perform one iterations of the Undecimated Transform and plot the result
% We use the daubechies series fo wavelets
[MLow3,MHigh3] = UWT(Spectra(:,2),Daub(2));
figure;
plot(Spectra(:,1),Spectra(:,2));
hold on
plot(Spectra(:,1),MHigh3);
title('2,4,6-trimethylphenol and one iteration of Daub(4)');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

% We first perform this peak detection on one iteration of the tranform

posneg = sign(MHigh3);
index = []; count = 1;
for i = 1:A-1
    if posneg(i) == 1 && posneg(i+1) == -1
        index(count) = i;
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
title('2,4,6-trimethylphenol and peak detection 2 iterations');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

%% Now perform peak detection algorithm on muliple iterations of Udecimated Wavelet Transform

% Perform 50 iterations of the Undecimated Transform and plot the result
[MLow4,MHigh4] = UWT1D(Spectra(:,2),Daub(2),50);

signdetect = [];
for i = 1:50
    posneg = sign(MHigh4(:,i));
    signdetect = [signdetect, posneg];
end

% This finds the peaks at each iteration of the transform
index = [];
for i = 50:-1:1
    count = 1;
    for j = 1:A-1
        if signdetect(j,i) == 1 && signdetect(j+1,i) == -1
            index(count,i) = j;
            count = count + 1;
        end
    end
    b(i) = count-1;
end

% Find the closest index values
Nindex = [];
count = 1;
for i = 1: b(50)
    col = 50;
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

% Plot the spectra and the found peaks
figure;
plot(Spectra(:,1),Spectra(:,2));
hold on
scatter(Npeaks(:,1),Npeaks(:,2));
title('2,4,6-trimethylphenol and peak detection 50 iterations');
set(gca,'xdir','reverse') % Flip the xscale to go from biggest to smallest

%% Conclusion: The Daubechies Series of Wavelets fail to smooth properly
        