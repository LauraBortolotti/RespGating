function [Trigger_sec, data, data_1, data_2] = RespGating_tq(MoPar, Time)
%
% Function description
% This function take the whole dataset of motion data.
% Dataset is: n-row as per time points acquired, 6-columns as per 6 motion
% parameters (Tx[m], Ty[m], Tz[m], Rx[deg], Ry[deg], Rz[deg])
%
% --------- Preprocessing:
% Convert T from [m] to [mm].
% live: data are acquired in [m]
%
% --------- First part:
    % I.a) It then select the best motion parameters to be used to track the respiration.
    % I.b) Then, it flip the motion parameters in order to have the exhale at the "bottom".
    % I.c) The number of time points corresponding to 1 resp cycle is found
%
% Suggestion for the live version: make it recognised by the acquisition of
% data at hold the breath condition.
%
% --------- Second part:
    %"Min" (end of the exhale) are identified
%
% --- Description of the variables used:
% Motion_Parameters = cell of strings reporting motion parameter names. They corresponds to the order to the 6-columns.
% MoPar = Motion parameters acquired using the motion camera. n-row as per time points acquired, 6-columns as per 6 motion.
% Trigger_sec = indexes at which found the gated data, n-row as per time points acquired, column vector
% data_MoPar = MoPar in [mm] and [deg]
% data = data_MoPar of the first 2 resp cycles (7 seconds)
% std_data = 6 std of data aka the deviation standard of each motion parameter.
% I = 6 indexes corresponding to the descend order of the std_data
% Motion_Parameters_sort = Motion_Parameters sorted by I
% data_norm = data normalised between 0 and 1
% mult = multiplicatin factor
% data_inv = data inverted by the multiplication factor
% data_inv_1 = data_inv is re-normalised between 0 and 1
% RespCycle = number of points corresponding to 1 resp cycle [number of points]
% RespCycle_points = 90% of the resp cycle in points [number of points]
%
% trigger = counter of the number of triggers found
% dt = average interval time between two consecutive data [seconds]
% RespCycle_sec = min lenght of the resp cycle from literature [seconds]
% Trigger_sec = time at wich the trigger is sent
% Trigger_Value = amplitutde of the predictor when trigger is sent
% flag_top = [0,1], it flag at wich phase of the resp data are. Goes to 1 when signal higher than running mean after last trigger (it goes over the peak)
% timesince = initial value of timesince last trigger
% istart_plot = [points]% initial data clean up/removal for plotting
%
% data = sum of the raw data, multiplied for the multiplicatin factor found
% data_1 = Perform a low pass filter using the previus 8 data
% data_2 = prospective drift removal
%
%
% --- Nested function (at the end of the document)
% AverageVsMedian
% ZerosVsOnes
% CharacteriseResp
%% Developping the code
%{
used during the developping:
load('C:\Users\ppzlb2\OneDrive - The University of Nottingham\UoN Upright scanner\Data\Body Tracking\LongDataset\Data_12-08-2022_A_long\Pose8.mat')
MoPar = [[DATA.RigidBodies_4.TX],...
    [DATA.RigidBodies_4.TY],...
    [DATA.RigidBodies_4.TZ],...
    DATA.RigidBodies_4.RX,...
    DATA.RigidBodies_4.RY,...
    DATA.RigidBodies_4.RZ];
Time = bsxfun(@minus,DATA.Time,DATA.Time(1,1));

RespGating(MoPar,Time)

%}
%
%% Useful variable
Motion_Parameters = {'Tx','Ty','Tz','qx','qy','qz','qw'};
%% Pre-processing

% Convert (Tx[m], Ty[m], Tz[m]) -NO-> (Tx[mm], Ty[mm], Tz[mm])
data_MoPar = [MoPar(:,1),...
    MoPar(:,2),...
    MoPar(:,3),...
    MoPar(:,4),...
    MoPar(:,5),...
    MoPar(:,6),...
    MoPar(:,7)];

% Live version: Have a plot in the GUI that updates while data are coming in.
% Obs: To take the changes in Motion parameters from the first timepoint lead to find the "min" relative to the 1st time point - don't do it.

% figure;plot(MoPar)

%% Part I
% From: C:\Users\ppzlb2\OneDrive - The University of Nottingham\UoN Upright scanner\Data\Body Tracking\BodyTracking_Exhale_ShortDataset_ISMRM2023_V2.m

% --- I.a)
% Consider the motion data recorded in the first 2 resp cycles of the data: find(Time>7,1)
% In the live version, this could be a "stop to acquire the data when 7 seconds are passed"
data = data_MoPar(2:find(Time>7,1),:);

% Deviation standard of the motion parameters is then found to rank them:
% std of data:
std_data = std(data,[],1);
% sort the std to rank them:
[~,I] = sort(std_data,'descend');
% I contains the index of the column

% The top 2 motion parameters showing the higher std are considered.
% Why 2? They will be summed-up. This will hilight the oscillatory
% behaviour ..but also the noise. We'll deal with the noise later in the
% code.

% "I(1:2)" select the first twos to be the parameters used for the prediction
% data, Motion_Parameters are then sorted based on it.
data = data(:,I(1:2)); % all the rows, the 2 columns with higher std 
Motion_Parameters_sort = Motion_Parameters(I);

% --- I.b)
%{
This part is to find out if motion data needs to be "flipped" (swap sign).

We expect data to show a sin (or cos) behaviour, having a "plateau" during exhale.
This influence the statistical caracterisation of the data:
- The average will be influnced on the numerical value of the data.
- The median on the occurrence of the numerical values.
So, by comparing them, we can find out the value of the plateau and the multiplication factor
(+1 if doesn't need to be inverted, -1 if it needs to be inverted).
This is done by using the nested function "AverageVsMedian".

Then, data series are sum to check if they are in phase or out of phase. As
data have been normalised:
- In case they are in pahse, the sum will be ~ 2
- In case they are out of phase, the sum will be ~ 0
The treshold is set as "<1.6" to identify if they are out of phase.

%}

% First data are normalised. This doesn't influence the relation with
% average and median and it is used to chek the in-phase/ou-of-phase later.
data_norm = normalize(data, 'range',[0 1]);

% Comparison between median and average
mult(1,1) = AverageVsMedian(data_norm(:,1));
mult(1,2) = AverageVsMedian(data_norm(:,2));

% Invert the data:
data_inv(:,1) = mult(1,1)*data_norm(:,1);
data_inv(:,2) = mult(1,2)*data_norm(:,2);

% Because the multiplication factor could be "-1", data could be now normalised 
% between [-1 0] instead of 0 and 1. So, it is re-normalise between 0 and 1
data_inv(:,1) = normalize(data_inv(:,1), 'range',[0 1]);
data_inv(:,2) = normalize(data_inv(:,2), 'range',[0 1]);


%{
Plot for debugging: in the live version, it could be an optional plot to
display.

fig = figure('Name','Data','Visible','on','Position',[0 0 1000 500]);

subplot(3,2,1)
hold on; plot(data_inv(:,1),'Displayname',['Rest' Motion_Parameters_sort{1,1}]);
yline(median(data_inv(:,1)),'Displayname',['mode' Motion_Parameters_sort{1,1}]);
yline(mean(data_inv(:,1)),':b','Displayname',['average' Motion_Parameters_sort{1,1}]);
legend;title([ ' : ' num2str(mult(1,1)) '*' Motion_Parameters_sort{1,1}])

subplot(3,2,3)
hold on; plot(data_inv(:,2),'Displayname',['Rest' Motion_Parameters_sort{1,2}]);
yline(median(data_inv(:,2)),'Displayname',['mode' Motion_Parameters_sort{1,2}]);
yline(mean(data_inv(:,2)),':b','Displayname',['average' Motion_Parameters_sort{1,2}]);
legend;title([ ' : ' num2str(mult(1,2)) '*' Motion_Parameters_sort{1,2}])
subplot(3,2,5)
hold on; plot(data_inv(:,1)+data_inv_1(:,2),'Displayname',['Rest' Motion_Parameters_sort{1,1} '+' Motion_Parameters_sort{1,2}]); 
legend; title('sum of the above')
%}


%{
If the data are in phase, the max of the sum of the motion
parameters is then ~ 2 as they are normalised between 0 and 1
if they are out of phase, means that the above code has not
working for one of them. Which one? Let's find out!

First of all, the external if function will compare the max of the sum with
the treshold value of 1.6 (the value 2 then is not used because it could 
set a too restrictive treshold. Data are noisy, so the max of two normalised
serie could be less than 2 as noise cancle out).

In case the sum is less than 1.6 then, the data are rounded.
As data are normalised, the effect of the function round() will be to bring
them either equal to 0 or to 1.

This will produce a squared data serie where "plateau" could be either
represented by the number of 0 (if it is at the bottom, where we'd like it
to be) or by a row of 1 (if data needs to be inverted again).

%}

if max(sum(data_inv,1))<1.6
    disp('Out of phase')
    
    % round() approximates the number to the nearest integer, 0 or
    % 1 as data are normalised.
    data_inv_norm_round(:,1) = round(data_inv(:,1));
    data_inv_norm_round(:,2) = round(data_inv(:,2));

    % As it has been done above:
    mult(2,1) = ZerosVsOnes(data_inv_norm_round(:,1));
    mult(2,2) = ZerosVsOnes(data_inv_norm_round(:,2));

    data_inv(:,2) = mult(2,2)*data_norm(:,2);
    data_inv(:,1) = mult(2,1)*data_norm(:,1);

    data_inv(:,1) = normalize(data_inv(:,1), 'range',[0 1]);
    data_inv(:,2) = normalize(data_inv(:,2), 'range',[0 1]);

    %{
    Plot for debugging: in the live version, it could be an optional plot to
    display.

    subplot(3,2,2)
    hold on; plot(data_inv_norm2(:,1),'Displayname',['Rest' Motion_Parameters_sort{1}]);
    plot(data_inv_norm_round(:,1),'Displayname',['round(Rest' Motion_Parameters_sort{1} ')']);
    legend;title([' : ' num2str(mult(1,1)) '*' Motion_Parameters_sort{1}])
    subplot(3,2,4)
    hold on; plot(data_inv_norm2(:,2),'Displayname',['Rest' Motion_Parameters_sort{2}]);
    plot(data_inv_norm_round(:,2),'Displayname',['round(Rest' Motion_Parameters_sort{2} ')']); 
    legend;title([' : ' num2str(mult(1,2)) '*' Motion_Parameters_sort{2}])
    subplot(3,2,6)
    hold on; plot(data_inv_norm2(:,1)+data_inv_norm2(:,2),'Displayname',['Rest' Motion_Parameters_sort{1} '+' Motion_Parameters_sort{2}]); 
    legend; title('sum of the above')
    %}

end

% figure;plot(Time(2:find(Time>7,1),:),data_inv)

%{
IN THE LIVE VERSION:
Here, the code has to stop, visualise the data in a plot and ask for the
approval for the scanner operator.
%}

% --- I.c)
% Evaluate the number of points coresponding to one resp cycle
RespCycle = CharacteriseResp(sum(data_inv,2),Time(2:find(Time>7,1)));
% Respiratory cycle varies during the period of the MRI scan.
% So, we are going to consider the 90% of the value found using the FFT
RespCycle_points = round(0.9*RespCycle); %rough min cycle length in POINTS

% clear variables no longer used
clear data data_inv data_inv_1 data_norm std_data Motion_Parameters


% In the live function: these will be new data coming in one by one and so
% Based on the results of part 1:
% Select only the motion parameters wanted
data = data_MoPar(:,I(1:2));
% Multiply for the multiplication factors in sequence as above
for m = 1:size(mult,1)
    data(:,1) = data(:,1).*mult(m,1); 
    data(:,2) = data(:,2).*mult(m,2);
end
clear m mult I
% Sum the signals:
data = sum(data,2);
% "data" is now the predictor to be used for finding the beginning of the exhale


%% Part II
% C:\Users\ppzlb2\OneDrive - The University of Nottingham\UoN Upright scanner\Data\Body Tracking\BodyTracking_Prospective_Longdataset_ISMRM2023.m

% Inizialise variables:
trigger = 0; % counting number of triggers found

% Average interval time between consecutive time points:
dt = mean(diff(Time));

% From the literature, we know that the average lenght of the resp cycle is
% 3 s, so we can consider 2 s as the minimum (this might need to be optimised)
RespCycle_sec = 2;

% Flag at wich moment of the resp cycle data is
flag_top = 0; % goes to 1 when signal higher than running mean after last trigger.. going over the peak

% Initial value for time since last trigger
timesince = 1000; % [sec]

% initial data clean up/removal for plotting
istart_plot = 4*RespCycle; % [points]

% Gradients: we expect to find the begin of the Exhale period when gradient
% is negative for at least a quarter of a cycle
grdlimL = fix(RespCycle_sec/4/dt); % [sec] distance over which grad must have been negative for a trigger (about a quarter of a cycle)

grdlimS = fix(grdlimL/4); % [sec] distance over which grad must have ticked up for a trigger @@4 needs optimizing



% Pre-allocate space for variables:
Trigger_sec = zeros(1,500); % Time at which triggers is sent
Trigger_Value = zeros(1,500); % Amplitude at triggers for plotting
% 500 is the number of triggers expected to be found, it should be set

data_2 = zeros(size(data,1),size(data,2));
data_1 = zeros(size(data,1),size(data,2)); 


% Live version: a number of data equal to the RespCycle needs to have been acquired
% before to start the anaylis. They could be the ones used in part 1: find(Time>7,1)

for j = find(Time>7,1)+1:length(Time) % j = RespCycle+1:length(Time)
    % Perform a low pass filter using the previus 8 data
    data_1(j) = sum(data(j-8:j))/9;
    % Peform a drift removal
    local_mean(j) = mean(data_1(j-RespCycle:j));data_2(j) = data_1(j)-local_mean(j);  %prospective drift removal
    %ALT h=fit(t(1:i),double(u1(1:i)),'poly3','Normalize','on'); %alternative method of prospective drift removal
    if j>istart_plot 
        running_min(j) = min(data_2(istart_plot:j)); running_mean(j) = mean(data_2(istart_plot:j)); running_max(j) = max(data_2(istart_plot:j));
    else %initial data clean up/removal for plotting- these parameters will be used to set axes
        running_min(j) = 0; running_mean(j) = 0; running_max(j) = 0;
    end
    local_min(j) = min(data_2(j-RespCycle:j)); local_mean(j)=mean(data_2(j-RespCycle:j)); local_max(j)=max(data_2(j-RespCycle:j));
    if (data_2(j)>local_mean(j)) % gone over peak- not reset to zero until a trigger point is discovered
        flag_top = 1; 
    end
    if (j<=grdlimL) grdlimL = RespCycle; end   
    if (j<=grdlimS) grdlimS = RespCycle; end   
    gradLong = data(j)-data(j-grdlimL); gradShort = data(j)-data(j-grdlimS);  % gradient must be negative and uptick recently for a trigger (with GradLong this determines early on flat part
    if trigger>0 %avoid triggers too soon, but also force inialization
        timesince = Time(j)-Trigger_sec(trigger);
    end
    if (flag_top==1 && gradLong<0 && gradShort>gradLong/15 && timesince>RespCycle_sec && data_2(j)<local_mean(j) && abs(data_2(j)-local_mean(j))>0.4*abs(local_min(j)-local_mean(j))  ) 
        trigger = trigger+1; %trigger count
        Trigger_sec(trigger) = Time(j);  Trigger_Value(trigger) = data_2(j);  
        flag_top = 0; %reset the flag the identifies going over top again before next trigger
    end
end
% %      axis([t(MinDel+1) t(temp) running_min(i) running_max(i)]); %taken out of loop as it slows the program.. but move it in loop for prospective? 
% name=strcat('Volunteer')
% Rabbit = figure('Name', name, 'Position', [500 0 1000 500])%, 'Visible','off');
% ax2 = subplot(2,1,2); hold(ax2,'on');
% ax1 = subplot(2,1,1); hold(ax1,'on');
% 
% % plot(ax1, Time,(x-x(1))/max(abs(x-x(1))));
% % plot(ax1, Time,(y-y(1))/max(abs(y-y(1)))); 
% % plot(ax1, Time,(z-z(1))/max(abs(z-z(1)))); 
% % plot(ax1, Time,(rx-rx(1))/max(abs(rx-rx(1)))); 
% % plot(ax1, Time,(ry-ry(1))/max(abs(ry-ry(1)))); 
% % plot(ax1, Time,(rz-rz(1))/max(abs(rz-rz(1)))); 
% plot(ax1, Time,(MoPar-MoPar(1,:))./max(abs((MoPar-MoPar(1,:))),[],2));
% 
% % plot(ax2, Trigger(find(Trigger~=0)),TriggerV(find(Trigger~=0)),'r*');
% % plot (ax2, Time,u2,'b'); 
% plot(ax2, Time,local_mean,'m'); 
% 
% % title([Volunteer{1,v} ' - ' Poses_Labelled{1,p} ' - Sternum']);
% linkaxes([ax1 ax2],'x')
% %xlim(ax1,[t(122) t(end)]); xlim(ax2,[t(122) t(end)-300]);
% axis(ax2,[Time(MinDel+1) Time(temp) running_min(i) running_max(i)]); %taken out of loop as it slows the program.. but move it in loop for prospective?        




%% Out of this programme:
% C:\Users\ppzlb2\OneDrive - The University of Nottingham\UoN Upright scanner\Data\Body Tracking\BodyTracking_Velocity_LongDataSet_ISMRM2023.m

end % end of the main function

%% Nested functions
% Obs: name of the variable have to be different from the one used in the
% main function because "global variable" screw things up in Matlab

function Multiplication_factor = AverageVsMedian(DATA_NORM)
% Find the multiplication factor by comparing the averange and the median.
% if the median is bigger than the average, it means that the number of
% occurrence of the plateau are happening near to the max instead of the
% min as we want
    if median(DATA_NORM(:,1)) > mean(DATA_NORM(:,1))
        Multiplication_factor = -1; % mode corrispond to the plateau
    else
       Multiplication_factor = +1;
    end
end

% --------------------------------------------------------------------------------------
function Multiplication_factor = ZerosVsOnes(DATA_INV_NORM_1_ROUND)
% Find the multiplication factor by comparing the number of zeros and ones.
% if number_of_zeros > number_of_ones -> the plateau is at 0,
% data don't need to be inverted
% if number_of_zeros < number_of_ones -> the plateau is at 1,
% data need to be inverted   

    if length(DATA_INV_NORM_1_ROUND(:,1))-nnz(DATA_INV_NORM_1_ROUND(:,1)) > nnz(DATA_INV_NORM_1_ROUND(:,1))
        % nnz() = number of non-zero element in an array
       Multiplication_factor = +1; % mode corrispond to the plateau
    else
       Multiplication_factor = -1;
    end
end

% -----------------------------------------------------------------------------------------
function [respCycle_frames] = CharacteriseResp(x,T,varargin)
% C:\Users\ppzlb2\OneDrive - The University of Nottingham\UoN Upright scanner\Data\Function\NUFFT\CharacteriseResp.m

% This function take data (x) and relative time string(T) and compute the
% nufft on apodised data (manipulate with exponential function). Time is
% elongated accordgly. It evaluate the main frequency found. It is
% optimised to find resp frequency for adult (0.1 < f < 0.4 Hz). If it is
% not found, then the default value is assigned (f = 0.3 [Hz]).

% Then, the relative lenght of the resp in seconds and number of frame is
% found. 

% Also, the trough and the  

%%
T = transpose(T);

    % (1) Create a longer time string uses to nufft(apodised data, long time string)
    t2 = [T, ...
    T+(T(end)-T(1)), ...
    T+ 2*(T(end)-T(1)), ...
    T+ 3*(T(end)-T(1))];
    
    X = bsxfun(@minus, x, mean(x,1));
    % (2) Apodization then zero filling
    x2 = zeros(size(t2,2),1); % Zero filling
    x2(1:size(X,1)) = X.*exp(-T.'/20); % Apodization
    % 1/20, 1/40, 1/80 will cut off the exponential ad 1/e, 1/e^2, 1/e^3
    clear X
    
    % Find the resp peak usinf nufft
    Y = nufft(x2,t2);
    Y = abs(Y);
    % Find the frequency range that is evaluated:
    nd = length(t2);
    f = (0:nd-1)/nd;
    % Find the resp frequency
    % For adults, respiration frequency is between 0.1 and 0.3 Hz
    Ymax = max(Y([find(f>=0.1,1,'first') : find(f>=0.4,1,'first')]));
    % Ymax might not be found, so:
    if isempty(Ymax) || Ymax == 0|| f(find(Y==Ymax))<0.2 || f(find(Y==Ymax)) >= 0.4
        disp('Max of the spectra not found OR smaller than 0.2 Hz OR greater than 0.4 Hz: f_resp will be the default one.')
        f_resp = 0.3; % default value is the average adult resp freq
    else
        f_resp = f(find(Y==Ymax));
    end
%     clear x2 t2
    
% Define the expected lenghts of one resp cycle:
    respCycle_seconds = 1/(f_resp); % [s]
    respCycle_frames = find(T>=(T(1,1)+respCycle_seconds),1,'first'); % [number of frames]                
       % NB: if T(1,1) is not 0, the above line might needs to be modified

    % Expected number of point withihn one respiration cycle
    % this number won't be costant over respiration cycle, and it
    % is used as an indication to find the following p-trough.
    % It is wise to underestimate it a bit, 80% of 2*HalfCycle
    

    if ~isempty(varargin)
        disp('Figure will be plotted')
        
        figure('Name','Apodisated data');
        plot(t2,x2,'Displayname','Apodised data');
        title('Apodisated data');
        xlabel('t [s]');
        ylabel('qw');

        figure('Name','Freq spectra');
        plot(f,Y,'Displayname','Freq Spectra');
        xlabel('f [Hz]'); ylabel('a.u');
        title('NUFFT (up to 1 Hz)'); legend;
        xline(f_resp,':k',{'Respiration'},'Displayname','f resp')

    end
end
