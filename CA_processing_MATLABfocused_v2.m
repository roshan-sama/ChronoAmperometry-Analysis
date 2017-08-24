% TO DO: Make Tau detection more robust?
% Subtract baseline that exists during prestep delay (before scaling) DONE
close all;
clear;
filesInDirectory = dir; %Struct containing properties of ALL files in directory
%Get base naming convention
baseName = input('Enter the base naming convention -this excludes the stage name \n(all the text before the step characters, for example, \nif the files were like C1_cortisol_ss_DSP.txt, then enter C1_cortisol_ss_\n','s');
fileNames_toBeProcessed = {}; %initialize file names cell. dir contains a list of ALL files, and needs to be curated to ensure that only data files are processed % The base naming convention will help in this process

%Get the names of the files that agree with the naming convenvtion in the directory
for k = 1:length(filesInDirectory)
    if strfind(filesInDirectory(k).name,baseName)           %Only files with the naming convention are selected     
        fileNames_toBeProcessed = [fileNames_toBeProcessed, filesInDirectory(k).name]; % variable remains a Cell each iteration
    end
end
%%
% Creating a new cell with strings containing only Stage and Trial name of each file
ExtractedStagePlusTrialNames = strrep(fileNames_toBeProcessed,baseName,''); 
ext = getFileExtension();
ExtractedStagePlusTrialNames = strrep(ExtractedStagePlusTrialNames,ext,'');
%FileProcessing variables: matrices
stageCurrentsMatrix         = [];   % Hold the currents corresponding to one stage
currentAveragesMatrix       = [];   % Place this Matrix in averages Sheet 
currentStdDevMatrix         = [];   % Place this Matrix in Std Dev Sheet
combinedData                = [];   % Place this Matrix in Combined Data Sheet. NOT SUPPORTED IN CURRENT VERSION
%Conglomerate variables: headers and names
StageNamesHeader_AvgStD     = {};   %Stage names header to be used in Average and Std Dev Matrix sheets
currentFileName             = {};
currentStageNamePlusTrial   = {};   % Cell for string value
currentStage                = {};   % Cell for string value
previousStage               = {};   % Cell for string
oneStageDoneFlag            = 0;    %logical

% Process files
for k = 1:length(fileNames_toBeProcessed)                           %Cycle through all valid files starting from first 
    currentFileName             =   fileNames_toBeProcessed(k);             %Obtain new file name      
    currentStageNamePlusTrial   =   ExtractedStagePlusTrialNames(k);
    if(contains(string(fileNames_toBeProcessed(k)),{'EIS','Ab','DSP','PBS'})) 
        continue;   %Ignore preliminary steps
    end        
    currentStage = extractStageName(currentStageNamePlusTrial);    
    newCurrent   = getArrayFromFile(string(currentFileName),'Current');   %Obtain new current data from file
       
    if(strcmp(currentStage,previousStage)) %If present file is not a new stage
        stageCurrentsArray = [stageCurrentsArray,newCurrent];% Augment current array with new currents        
        previousStage = currentStage;  %The previous stage used in the next run is set here
    else %New Stage processing             
        if(oneStageDoneFlag) %Only average the previous stage if there exists a previous stage!            
            currentAveragesMatrix   =   [currentAveragesMatrix mean(stageCurrentsArray,2)];            
            currentStdDevMatrix     =   [currentStdDevMatrix std(stageCurrentsArray,0,2)];        
        else
            oneStageDoneFlag = 1;
        end
        StageNamesHeader_AvgStD =   [StageNamesHeader_AvgStD currentStage];       %Add the new stage's name to the header 
        previousStage = currentStage;       % Set previous stage for new run
        stageCurrentsArray = newCurrent;    % Refresh current array and place newCurrent array in it
    end      
end

currentAveragesMatrix   =   [currentAveragesMatrix mean(stageCurrentsArray,2)]; %Process the last stage's currents   
currentStdDevMatrix     =   [currentStdDevMatrix std(stageCurrentsArray,0,2)];
timeArray = getArrayFromFile(string(currentFileName), 'Time' );
voltageArray = getArrayFromFile(string(currentFileName),'Voltage');

%Data Entry 
%Reordering Columns if needed %Cnt DosesArray = [0,1,2,3,4,6,8,10,15,20,25,30,35,40,50,70,80,90,100,120,135,150];
%For continuous dosing: reorderArray = [1,19,8,20,21,22,9,2,10,11,3,12,13,4,14,5,15,6,16,17,7,18];
%For BEV: reorderArray = [3,2,1];
% reorderArray = [1,19,8,20,21,22,9,2,10,11,3,12,13,4,14,5,15,6,16,17,7,18];
% StageNamesHeader_AvgStD = reorderColumns( StageNamesHeader_AvgStD, reorderArray );
% currentAveragesMatrix = reorderColumns( currentAveragesMatrix, reorderArray );
% currentStdDevMatrix = reorderColumns( currentStdDevMatrix, reorderArray );
%Writing the data to file
xlswrite('Data Analysis File',StageNamesHeader_AvgStD,  'Current Averages','A1')
xlswrite('Data Analysis File',currentAveragesMatrix,    'Current Averages','A2')
xlswrite('Data Analysis File',StageNamesHeader_AvgStD,  'Std Devs','A1')
xlswrite('Data Analysis File',currentStdDevMatrix,      'Std Devs','A2')
xlswrite('Data Analysis File',timeArray,      'Combined Data','A2')

%Data Analysis and processing:
%Get the times at which steps occur
v_step1 = 0.27;v_step2 = -0.27;t_step1 = 0; t_step2 = 0; index_step1 = 0; index_step2 = 0;
for k = 1:length(timeArray) 
    if voltageArray(k) > v_step1 %If voltage is found to exceed 0.2 volts, then this is time of first step        
        t_step1 = timeArray(k);
        index_step1 = k;
        break;
    end
end
for k = 1:length(timeArray) 
    if voltageArray(k) < v_step2 %If voltage is found to exceed 0.2 volts, then this is time of first step        
        t_step2 = timeArray(k);
        index_step2 = k;
        break;
    end
end
%Variables to help and hold processed data
[r, c] = size(currentAveragesMatrix); %Rows and columns of the data Matrix
Tau = zeros(2,c);
scaledCurrents= zeros(r,c); %Initialize a scaled version of the currents matrix
for k = 1:c %The type of scaling here is as follows: the baseline is subtracted and the peak current is treated as unity, and all other values are ratio w.r.t it
    %Variables needed for data curation: Baseline, scaler
    baseline = currentAveragesMatrix(index_step1-1,k);
    scaledCurrents(:,k) = currentAveragesMatrix(:,k) - baseline; %Obtainin baseline adjusted values for current stage
    scaler1             = scaledCurrents(index_step1,k); %get the scaling factor for the first step of current stage     
    scaler2             = scaledCurrents(index_step2,k); %get the scaling factor for the second step of current stage
    scaledCurrents(index_step1-1:index_step2-1,k) = -1*scaledCurrents(index_step1-1:index_step2-1,k)/scaler; %apply scaling, stage 1
    scaledCurrents(index_step2:end,k) = -1*scaledCurrents(index_step2:end,k)/scaler; %apply scaling, stage 2
    Tau(:,k)            = getTimeConstant(timeArray,currentAveragesMatrix(:,k),t_step1,t_step2); %deduce Tau
end
xlswrite('Data Analysis File',StageNamesHeader_AvgStD,  'Tau','A1')
xlswrite('Data Analysis File',Tau,      'Tau','A2')
plot(timeArray,scaledCurrents)
legend(StageNamesHeader_AvgStD)
saveas(gcf,'Scaled Currents.fig')
figure;
plot(1:c,Tau)
saveas(gcf,'Tau.fig')


%-- END PROGRAM --%

%% Get the data file's extension from the user
function [fExt] = getFileExtension()%Get the file extension
    fExt = input('Enter the file extension (.dat, .txt, .xls .xlsx','s');    
    while ~mean(strcmp(fExt,{'.txt','.dat','.xls','.xlsx'}))
        fExt = input('Valid file Extensions are: .dat, .txt, .xls etc','s');
    end
end
%% Return a current or Time Array from ChronoAmp Data file
function [ ReturnArray ] = getArrayFromFile( fileName, returnArray )
%getTimeAndCurrentArrays- Returns time and current arrays of chronoamp data
%   files may be in  .DTA or .txt format
    id = fopen(fileName);
    a = textscan(id, '%s %s %s %s %s %s %s %s %s');
    ReturnArray = [];        
    startRow    = 1;    
    for k = 1:length(a{1,2}) %Find start row
        if strcmp(a{1,2}{k},'bits')
            startRow = k+1;
        end        
    end    
    if strcmp(returnArray,'Time') % Return either time, current or voltage
        lalalalalalasmurfAhappySonng = 2;
    elseif strcmp(returnArray,'Current')
        lalalalalalasmurfAhappySonng = 4;
    elseif strcmp(returnArray,'Voltage')
        lalalalalalasmurfAhappySonng = 3;
    end    
    for k = startRow:length(a{1,2}) % Collect entries into an array (column array)
        ReturnArray = [ReturnArray ; str2double(a{1,lalalalalalasmurfAhappySonng}{k})];
    end
    fclose(id); % Close the file
end
%% Extract stage name
function [ stageName] = extractStageName(stagePlusTrial)    
    Temp = char(stagePlusTrial);    %Convert and store the input cell string in a char array
    stageName = Temp(1:end-3);      %Remove the last 3 letters: _#1, _#2, etc.
end
%% Initialize variable to hold the index of time corresponding to the first potential step%%
function [ Tau ] = getTimeConstant( t, c, t_step1, t_step2 )   
    if(~isrow(t)) 
        time      = t';
    else
        time = t;
    end
    if(~isrow(c)) 
        current   = c';
    else
        current = c;
    end
    % Gathering needed variables ---------
    t_start_index1 = 1; %Initialize variable to hold the index of time corresponding to the first potential step    
    for t = time
        if(t >= t_step1)
            break;
        else
            t_start_index1 = t_start_index1 + 1;
        end
    end   % Obtain index of arrays coresponding to the first step
    I_start1 = current(t_start_index1); % Store Initial current corresponding to the first potential step    
    
    t_start_index2 = 1; %Initialize variable to hold the index of time corresponding to the second potential step    
    for t = time
        if(t >= t_step2)
            break;
        else
            t_start_index2 = t_start_index2 + 1;
        end
    end   % Obtain index or array coresponding to the SECOND step
    I_start2 = current(t_start_index2); % Store Initial current corresponding to the SECOND potential step    
    
    t_end_index1 = t_start_index2 - 1;   % Indices of the respective ends of each step
    t_end_index2 = length(time);
    
    % Detect minimum, this will be I final1    
    I_min1  = I_start1;
    for I = current(t_start_index1:t_end_index1)
        if(I < I_min1)
            I_min1 = I;        
        end
    end
    % Detect minimum, this will be I final2
    I_min2  = I_start2;
    for I = current(t_start_index2:t_end_index2)
        if(I > I_min2)
            I_min2 = I;        
        end
    end
    
    DelI_step1 = I_start1 - I_min1;
    DelI_step2 = I_start2 - I_min2;    
    
    I_target_Index1  = t_start_index1; %Obtain target index 1
    for I = current(t_start_index1:t_end_index1)
        if  I <= 0.36*DelI_step1 + I_min1
            break;
        else
            I_target_Index1 = I_target_Index1 + 1;
        end
    end   
    I_target_Index2  = t_start_index2; %obtain target index 2
    for I = current(t_start_index2:t_end_index2)
        if  I >= 0.36*DelI_step2 + I_min2
            break;
        else
            I_target_Index2 = I_target_Index2 + 1;
        end
    end  
    
    Tau1 = time(I_target_Index1)-time(t_start_index1); % Time difference obtained by taking difference between preceeding entities
    Tau2 = time(I_target_Index2)-time(t_start_index2); % Time difference obtained by taking difference between preceeding entities
    Tau = [Tau1;Tau2];
end
%% Function reorders a matrix's columns based on an input orderArray
function [ returnMatrix ] = reorderColumns( inputMatrix, orderArray )%Function reorders a matrix's columns based on an input orderArray
if iscell(inputMatrix)
    returnMatrix = {};  
    for k = 1:length(orderArray)
        returnMatrix(orderArray(k)) = inputMatrix(k);
    end
else
    returnMatrix = zeros(size(inputMatrix));

    for k = 1:length(orderArray)
        returnMatrix(:,orderArray(k)) = inputMatrix(:,k);
    end
end
end