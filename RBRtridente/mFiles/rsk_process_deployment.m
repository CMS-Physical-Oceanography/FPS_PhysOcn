function [RSK_data] = rsk_process_deployment(filename,RSKtoolspath)

%=========================================================================
% rsk_proccess_deployment reads the .rsk file output from the Ruskin.exe
% software. The inputs to the function are
%     - filename: the full filename ".rsk" 
%     - RSKtools: file directory to the RSKtools
%
%|RSKtools| is RBR's open source Matlab toolbox for reading,
% visualizing, and post-processing RBR logger data. It provides
% high-speed access to large RSK data files. Users may plot data as a
% time series or as depth profiles using tailored plotting
% utilities. If the tool box is not already downloaded, it
%can be found at this link: https://rbr-global.com/support/matlab-tools/
%
% The output for this file is a data structure that contains time,
% chlorophyll-a, FDOM, and turbidity data and units. 
%
%If a user needs more clarification on certain parts of the tool box type
%"help RSKtools" in the command window. 
%=========================================================================

%Path to tool box RSKtools
addpath(RSKtoolspath);

rsk=RSKopen(filename);

startTime=rsk.epochs.startTime;
endTime=rsk.epochs.endTime;

% Specify the time frame to read data from the RBR
RSK=RSKreaddata(rsk,'t2',startTime+1);
%rsk = RSKreaddata(rsk, 't1', t1, 't2', t2);

% RSKplotdata(RSK)

%Data Output Structure
RSK_data.chlorophyll_a.data=RSK.data.values(:,1);
RSK_data.chlorophyll_a.units='ug/l';
RSK_data.FDOM.data=RSK.data.values(:,2);
RSK_data.FDOM.units='ppb';
RSK_data.turbidity.data=RSK.data.values(:,3);
RSK_data.turbidity.units='FTU';
RSK_data.time=RSK.data.tstamp;

figure(); clf
subplot(311)
    plot(RSK_data.time, RSK_data.chlorophyll_a.data); grid on; 
    set(gca,'XTickLabel',[]); ylabel(RSK_data.chlorophyll_a.units); 
    title('Raw RBR Tridente Data');
subplot(312)
    plot(RSK_data.time, RSK_data.FDOM.data); grid on; 
    set(gca,'XTickLabel',[]); ylabel(RSK_data.FDOM.units);
subplot(313)
    plot(RSK_data.time, RSK_data.turbidity.data); grid on; 
    ylabel(RSK_data.turbidity.units);









