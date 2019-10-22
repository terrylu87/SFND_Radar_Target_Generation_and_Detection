clear all
clc;
 
%% Radar Specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%DONE* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
target_position = 60; %distance from the sensor
target_velocity = 12; % m/s
 
 
%% FMCW Waveform Generation
 
% *%DONE* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.
 
%Operating carrier frequency of Radar
fc= 77e9;              %carrier freq
max_range = 200;       %200m
range_resolution = 1;  %1m
max_velocity = 70;     %70m/s
c = 3e8;               %speed of light
Bsweep = c / (2*range_resolution); %bandwidth
Tchirp = 5.5*2*max_range/c;
Slope = Bsweep / Tchirp;
                                                         
%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation.
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps
 
%The number of samples on each chirp.
Nr=1024;                  %for length of time OR # of range cells
 
% Timestamp for running the displacement scenario for every sample on each
% chirp
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples
 
 
%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal
 
%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));
 
 
%% Signal generation and Moving Target simulation
% Running the radar scenario over the time.
 
for i=1:length(t)        
   
   
    % *%DONE* :
    %For each time stamp update the Range of the Target for constant velocity.
    r_t(i) = target_position + target_velocity * t(i);
    td(i) = 2 * r_t(i) / c;
    % *%DONE* :
    %For each time sample we need update the transmitted and
    %received signal.
    Tx(i) = cos(2*pi * ( fc*t(i) + Slope * t(i)^2 / 2));
    Rx(i) = cos(2*pi * ( fc*(t(i)-td(i)) + Slope * (t(i)-td(i))^2 / 2));
   
    % *%DONE* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i)*Rx(i);
   
end
 
%% RANGE MEASUREMENT
 
 
 % *%DONE* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix,[Nr,Nd]);
 % *%DONE* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
signal_fft = fft(Mix,Nr);
signal_fft = signal_fft./Nr;

 % *%DONE* :
% Take the absolute value of FFT output
signal_fft = abs(signal_fft);
 % *%DONE* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
signal_fft = signal_fft(1:Nr/2); 
 
%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)
 
 % *%DONE* :
 % plot FFT output
plot(signal_fft);
 
axis ([0 200 0 1]);
 
 
 
%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM
 
 
% Range Doppler Map Generation.
 
% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.
 
Mix=reshape(Mix,[Nr,Nd]);
 
% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);
 
% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;
 
%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);
 
%% CFAR implementation
 
%Slide Window through the complete Range Doppler Map
 
% *%DONE* :
%Select the number of Training Cells in both the dimensions.
Tr = 10;
Td = 8;

% *%DONE* :
%Select the number of Guard Cells in both dimensions around the Cell under
%test (CUT) for accurate estimation
Gr = 4;
Gd = 4;

% *%DONE* :
% offset the threshold by SNR value in dB
offset = 8; %db

% *%DONE* :
% The threshold of each cell is the averange noise value from the training
% cells. Use a convolution with a customized kernel to produce the result.
% Here is an example of kernel which Tr=Gr=Td=Gd=1
% kernel:
% 1/16, 1/16, 1/16, 1/16, 1/16,
% 1/16, 0,    0,    0,    1/16,
% 1/16, 0,    0,    0,    1/16,
% 1/16, 0,    0,    0,    1/16,
% 1/16, 1/16, 1/16, 1/16, 1/16,

kernel = ones((Tr+Gr)*2+1, (Td+Gd)*2+1);
for i=(1+Tr):(1+Tr+Gr*2)
    for j=(1+Td):(1+Td+Gd*2)
        kernel(i,j)=0;
    end
end
kernel = kernel/sum(sum(kernel));
DynamicThreshold = pow2db(conv2(db2pow(RDM),kernel,'same')) + offset;

% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
% CFAR
RDM = double(RDM >= DynamicThreshold);

% *%DONE* :
% The process above will generate a thresholded block, which is smaller
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0.
RDM(union(1:(Tr+Gr),end-(Tr+Gr-1):end),:) = 0;  % Rows
RDM(:,union(1:(Td+Gd),end-(Td+Gd-1):end)) = 0;  % Columns 

% *%DONE* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
%figure,surf(doppler_axis,range_axis,'replace this with output');
figure,surf(doppler_axis,range_axis,RDM);
colorbar;
