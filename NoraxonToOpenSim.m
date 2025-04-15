% 3.15.2025
% Vicki Wang
% Human Engineering Research Laboratories

% This code converts Noraxon files into:
% 1 - marker data for the entire trial (.trc)
% 2 - marker data for a static pose (first 100pts (0.2s) of trial) (.trc)
% 3 - GRF and CoP data for the trial (.mot)

%%
clear
clc
close all

%% READ DATA FROM NORAXON CSV

% Change parameters here
subjectName = 'Subject1';
bodyweight = 180; % in POUNDS
trialCSV = 'ExampleNoraxonData.csv';
%

trialID = trialCSV(18:end-4);

% Change these as needed to pull files from your system
data_read = readtable(trialCSV);
load("/MATLAB Drive/NoraxonTrajectoriesAndForces.mat") % Loads trajectory and force names
trial_trc = readcell("marker_trcBLANK.txt");
trial_mot = readcell("grf_motBLANK.txt");
%

function trial_data = readNoraxon(trial_read, trajectories, forces)
    columnPrefix = "NoraxonMyoMotion_Trajectories_";
    trial_data.time = table2array(trial_read(:,'time'));
    for i = 1:length(trajectories)
        if contains(trajectories{i}, "Metatarso")
            columnSuffix = "_m";
        else
            columnSuffix = "_mm_";    
        end
        columnName = append(columnPrefix, trajectories{i}, columnSuffix);
        if contains(trajectories{i}, "1st") % Need to change first letter to a non-number otherwise MATLAB is unhappy
            trajectories{i}(1:3) = 'Med';
            trial_data.(trajectories{i}) = table2array(trial_read(:, columnName));
        elseif contains(trajectories{i}, "5th")
            trajectories{i}(1:3) = 'Lat';
            trial_data.(trajectories{i}) = table2array(trial_read(:, columnName));
        elseif contains(trajectories{i}, "7th")
            trajectories{i}(1:3) = 'Svn';
            trial_data.(trajectories{i}) = table2array(trial_read(:, columnName));
        elseif contains(trajectories{i}, "10th")
            trajectories{i}(1:3) = 'Ten';
            trial_data.(trajectories{i}) = table2array(trial_read(:, columnName));
        else
            trial_data.(trajectories{i}) = table2array(trial_read(:, columnName));
        end
    end

    for i = 1:length(forces)
        trial_data.(forces{i}) = table2array(trial_read(:, forces(i)));
    end
end

trial_data = readNoraxon(data_read, trajectories, forces);

%% FIX TRAJECTORY LABELS

% Need to change first letter to a non-number otherwise MATLAB is unhappy
for i = 1:length(trajectories)
    if contains(trajectories{i}, "1st") 
        trajectories{i}(1:3) = 'Med';
    elseif contains(trajectories{i}, "5th")
        trajectories{i}(1:3) = 'Lat';
    elseif contains(trajectories{i}, "7th")
        trajectories{i}(1:3) = 'Svn';
    elseif contains(trajectories{i}, "10th")
        trajectories{i}(1:3) = 'Ten';
    end
end

%% ROTATE DATA

% Switch y and z, then negate the new z
function rotated_data = rotateAroundX(trial_data, trajectories)
    rotated_data = trial_data;
    for i = 1:length(trajectories)
        if contains(trajectories{i}, '_y')
            markerName = erase(trajectories{i}, '_y');
            YZmarkers = {};
            YZmarkers{1} = append(markerName, '_y');
            YZmarkers{2} = append(markerName, '_z');
            tempY = trial_data.(YZmarkers{1});
            tempZ = trial_data.(YZmarkers{2});
            rotated_data.(YZmarkers{1}) = tempZ;
            rotated_data.(YZmarkers{2}) = -tempY;
        end
    end
end

rotated_data = rotateAroundX(trial_data, trajectories);

%% ROTATE AND TRANSLATE CENTER OF PRESSURE

% Noraxon insoles offset on the x axis by 1.8" - correct for this
% Also, change NaNs to zeros, and make corresponding GRF also zero
function [fixCoPx, fixCoPy, fixGRF] = fixCoPNaNs(rawCoPx, rawCoPy, rawGRF)
    fixCoPx = rawCoPx;
    fixCoPy = rawCoPy;
    fixGRF = rawGRF;
    for t = 1:length(rawCoPx)
        if isnan(rawCoPx(t))
            fixCoPx(t) = 0;
            fixCoPy(t) = 0;
            fixGRF(t) = 0;
        else
            fixCoPx(t) = rawCoPx(t) - 50;
            fixCoPy(t) = rawCoPy(t) - 50;
        end
    end
end

[rotated_data.noNaN_CoP_LT_x, rotated_data.noNaN_CoP_LT_y, rotated_data.noNaN_GRF_LT] = fixCoPNaNs(rotated_data.Insole_LTInsole_CenterOfPressure_x_mm_, rotated_data.Insole_LTInsole_CenterOfPressure_y_mm_, rotated_data.LTInsole_Total___);
[rotated_data.noNaN_CoP_RT_x, rotated_data.noNaN_CoP_RT_y, rotated_data.noNaN_GRF_RT] = fixCoPNaNs(rotated_data.Insole_RTInsole_CenterOfPressure_x_mm_, rotated_data.Insole_RTInsole_CenterOfPressure_y_mm_, rotated_data.RTInsole_Total___);

% Midpoint - x
footVector.frontLT_x = (rotated_data.MedMetatarsophalangealJointLT_x + rotated_data.LatMetatarsophalangealJointLT_x)/2;
footVector.frontLT_y = (rotated_data.MedMetatarsophalangealJointLT_y + rotated_data.LatMetatarsophalangealJointLT_y)/2;
footVector.frontLT_z = (rotated_data.MedMetatarsophalangealJointLT_z + rotated_data.LatMetatarsophalangealJointLT_z)/2;
footVector.frontRT_x = (rotated_data.MedMetatarsophalangealJointRT_x + rotated_data.LatMetatarsophalangealJointRT_x)/2;
footVector.frontRT_y = (rotated_data.MedMetatarsophalangealJointRT_y + rotated_data.LatMetatarsophalangealJointRT_y)/2;
footVector.frontRT_z = (rotated_data.MedMetatarsophalangealJointRT_z + rotated_data.LatMetatarsophalangealJointRT_z)/2;

% Center trajectories with the heel as the origin
footVector.frontLT_x_centered = footVector.frontLT_x - rotated_data.HeelBackLT_x;
footVector.frontLT_y_centered = footVector.frontLT_y - rotated_data.HeelBackLT_y;
footVector.frontLT_z_centered = footVector.frontLT_z - rotated_data.HeelBackLT_z;
footVector.frontRT_x_centered = footVector.frontRT_x - rotated_data.HeelBackRT_x;
footVector.frontRT_y_centered = footVector.frontRT_y - rotated_data.HeelBackRT_y;
footVector.frontRT_z_centered = footVector.frontRT_z - rotated_data.HeelBackRT_z;

% Find rotation matrices through each time point
for i = 1:length(rotated_data.time)
    footVecLT = [footVector.frontLT_x_centered(i); footVector.frontLT_y_centered(i); footVector.frontLT_z_centered(i)];
    footVecRT = [footVector.frontRT_x_centered(i); footVector.frontRT_y_centered(i); footVector.frontRT_z_centered(i)];
    CoPVecLT = [rotated_data.noNaN_CoP_LT_y(i); 0; rotated_data.noNaN_CoP_LT_x(i)];
    CoPVecRT = [rotated_data.noNaN_CoP_RT_y(i); 0; rotated_data.noNaN_CoP_RT_x(i)];
    normFootLT = norm(footVecLT);
    normFootRT = norm(footVecRT);
    footUnitVecLT = footVecLT/normFootLT;
    footUnitVecRT = footVecRT/normFootRT;
    insoleVec = [1; 0; 0];
    rotMatrixLT = footUnitVecLT/insoleVec;
    rotMatrixRT = footUnitVecRT/insoleVec;
    rotatedCoPLT = rotMatrixLT*CoPVecLT;
    rotatedCoPRT = rotMatrixRT*CoPVecRT;
    rotated_data.CoP_LT_x_centered(i) = rotatedCoPLT(1,1);
    rotated_data.CoP_LT_y_centered(i) = rotatedCoPLT(2,1);
    rotated_data.CoP_LT_z_centered(i) = rotatedCoPLT(3,1);
    rotated_data.CoP_RT_x_centered(i) = rotatedCoPRT(1,1);
    rotated_data.CoP_RT_y_centered(i) = rotatedCoPRT(2,1);
    rotated_data.CoP_RT_z_centered(i) = rotatedCoPRT(3,1);
end

% % Plot foot as it rotates around the heel, with rotated CoP mapped
% figure(1)
% for t = 1:50:length(rotated_data.time)
%     foot = [footVector.frontLT_x_centered(t), footVector.frontLT_y_centered(t), footVector.frontLT_z_centered(t); 0, 0, 0];
%     CoP = [rotated_data.CoP_LT_x_centered(t), rotated_data.CoP_LT_y_centered(t), rotated_data.CoP_LT_z_centered(t)];
%     toesMed = [rotated_data.FootToeLT_x(t) - rotated_data.HeelBackLT_x(t), rotated_data.FootToeLT_y(t) - rotated_data.HeelBackLT_y(t), rotated_data.FootToeLT_z(t) - rotated_data.HeelBackLT_z(t); rotated_data.MedMetatarsophalangealJointLT_x(t) - rotated_data.HeelBackLT_x(t), rotated_data.MedMetatarsophalangealJointLT_y(t) - rotated_data.HeelBackLT_y(t), rotated_data.MedMetatarsophalangealJointLT_z(t) - rotated_data.HeelBackLT_z(t)];
%     toesLat = [rotated_data.FootToeLT_x(t) - rotated_data.HeelBackLT_x(t), rotated_data.FootToeLT_y(t) - rotated_data.HeelBackLT_y(t), rotated_data.FootToeLT_z(t) - rotated_data.HeelBackLT_z(t); rotated_data.LatMetatarsophalangealJointLT_x(t) - rotated_data.HeelBackLT_x(t), rotated_data.LatMetatarsophalangealJointLT_y(t) - rotated_data.HeelBackLT_y(t), rotated_data.LatMetatarsophalangealJointLT_z(t) - rotated_data.HeelBackLT_z(t)];
%     plot3(foot(:,1), foot(:,2), foot(:,3), 'b', ...
%         CoP(:,1), CoP(:,2), CoP(:,3), 'ro', ...
%         toesMed(:,1), toesMed(:,2), toesMed(:,3), 'b', ...
%         toesLat(:,1), toesLat(:,2), toesLat(:,3), 'b')
%     view(-180,0)
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     xlim([-20 300])
%     ylim([-160 160])
%     zlim([-160 160])
%     title('Centered COP mapped over centered foot trajectories')
%     pause(0.05)
% end

% Translate CoPs to move with the foot through space
for t = 1:length(rotated_data.time)
    rotated_data.CoP_LT_x(t) = rotated_data.CoP_LT_x_centered(t) + rotated_data.HeelBackLT_x(t);
    rotated_data.CoP_LT_y(t) = rotated_data.CoP_LT_y_centered(t) + rotated_data.HeelBackLT_y(t);
    rotated_data.CoP_LT_z(t) = rotated_data.CoP_LT_z_centered(t) + rotated_data.HeelBackLT_z(t);
    rotated_data.CoP_RT_x(t) = rotated_data.CoP_RT_x_centered(t) + rotated_data.HeelBackRT_x(t);
    rotated_data.CoP_RT_y(t) = rotated_data.CoP_RT_y_centered(t) + rotated_data.HeelBackRT_y(t);
    rotated_data.CoP_RT_z(t) = rotated_data.CoP_RT_z_centered(t) + rotated_data.HeelBackRT_z(t);
end

% % Plot final CoP over original trajectory data
% figure(2)
% for t = 1:50:length(rotated_data.time)
%     footL = [footVector.frontLT_x(t), footVector.frontLT_y(t), footVector.frontLT_z(t); rotated_data.HeelBackLT_x(t), rotated_data.HeelBackLT_y(t), rotated_data.HeelBackLT_z(t)];
%     toesMedL = [rotated_data.FootToeLT_x(t), rotated_data.FootToeLT_y(t), rotated_data.FootToeLT_z(t); rotated_data.MedMetatarsophalangealJointLT_x(t), rotated_data.MedMetatarsophalangealJointLT_y(t), rotated_data.MedMetatarsophalangealJointLT_z(t)];
%     toesLatL = [rotated_data.FootToeLT_x(t), rotated_data.FootToeLT_y(t), rotated_data.FootToeLT_z(t); rotated_data.LatMetatarsophalangealJointLT_x(t), rotated_data.LatMetatarsophalangealJointLT_y(t), rotated_data.LatMetatarsophalangealJointLT_z(t)];
%     CoPL = [rotated_data.CoP_LT_x(t), rotated_data.CoP_LT_y(t), rotated_data.CoP_LT_z(t)];
%     %
%     footR = [footVector.frontRT_x(t), footVector.frontRT_y(t), footVector.frontRT_z(t); rotated_data.HeelBackRT_x(t), rotated_data.HeelBackRT_y(t), rotated_data.HeelBackRT_z(t)];
%     toesMedR = [rotated_data.FootToeRT_x(t), rotated_data.FootToeRT_y(t), rotated_data.FootToeRT_z(t); rotated_data.MedMetatarsophalangealJointRT_x(t), rotated_data.MedMetatarsophalangealJointRT_y(t), rotated_data.MedMetatarsophalangealJointRT_z(t)];
%     toesLatR = [rotated_data.FootToeRT_x(t), rotated_data.FootToeRT_y(t), rotated_data.FootToeRT_z(t); rotated_data.LatMetatarsophalangealJointRT_x(t), rotated_data.LatMetatarsophalangealJointRT_y(t), rotated_data.LatMetatarsophalangealJointRT_z(t)];
%     CoPR = [rotated_data.CoP_RT_x(t), rotated_data.CoP_RT_y(t), rotated_data.CoP_RT_z(t)];
%     %
%     plot3(footR(:,1), footR(:,2), footR(:,3), 'b', ...
%         toesMedR(:,1), toesMedR(:,2), toesMedR(:,3), 'b', ...
%         toesLatR(:,1), toesLatR(:,2), toesLatR(:,3), 'b', ...
%         CoPR(:,1), CoPR(:,2), CoPR(:,3), 'ro', ...
%         footL(:,1), footL(:,2), footL(:,3), 'r', ...
%         toesMedL(:,1), toesMedL(:,2), toesMedL(:,3), 'r', ...
%         toesLatL(:,1), toesLatL(:,2), toesLatL(:,3), 'r', ...
%         CoPL(:,1), CoPL(:,2), CoPL(:,3), 'bo')
%     %
%     view(-180,0)
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     xlim([-100 400])
%     zlim([-400 400])
%     title('CoP mapped over raw foot trajectories')
%     pause(0.05)
% end

% Convert from mm to m
rotated_data.CoP_LT_x = rotated_data.CoP_LT_x / 1000;
rotated_data.CoP_LT_y = rotated_data.CoP_LT_y / 1000;
rotated_data.CoP_LT_z = rotated_data.CoP_LT_z / 1000; 
rotated_data.CoP_RT_x = rotated_data.CoP_RT_x / 1000;
rotated_data.CoP_RT_y = rotated_data.CoP_RT_y / 1000;
rotated_data.CoP_RT_z = rotated_data.CoP_RT_z / 1000; 

%% CONVERT GRFS TO NEWTONS
rotated_data.GRF_LT = (rotated_data.noNaN_GRF_LT ./ 100) .*bodyweight.*4.448;
rotated_data.GRF_RT = (rotated_data.noNaN_GRF_RT ./ 100) .*bodyweight.*4.448;

%% CALCULATE GROUND REACTION MOMENTS
rotated_data.M_LT_x = rotated_data.GRF_LT.*rotated_data.CoP_LT_z;
% rotated_data.M_LT_y = 0;
rotated_data.M_LT_z = rotated_data.GRF_LT.*rotated_data.CoP_LT_x;
rotated_data.M_RT_x = rotated_data.GRF_RT.*rotated_data.CoP_RT_z;
% rotated_data.M_RT_y = 0;
rotated_data.M_RT_z = rotated_data.GRF_RT.*rotated_data.CoP_RT_x;


%% SET UP .TRC FILES

% Set up header info
trial_trcName = append(subjectName, '_', trialID, '.trc');
dataCount = length(trial_data.time);

%%

% Add specific header info
trial_trc{1,4} = trial_trcName;
trial_trc{3,3} = dataCount;
trial_trc{3,8} = dataCount;

% Add frame, time, and marker coordinate data
for i = 1:dataCount
    trial_trc{i+6,1} = i; % Frame#
    trial_trc{i+6,2} = i*0.002-0.002; % Time
    for k = 1:length(trajectories) % Marker coordinate data
        trial_trc{i+6,k+2} = rotated_data.(trajectories{k})(i);
    end
end

% Remove <missing> in cells
mask = cellfun(@(x) any(isa(x,'missing')), trial_trc);
trial_trc(mask) = {[]}; % or whatever value you want to use

%% WRITE MARKER TRC FILES

% Trial
trial_fileName = trial_trcName(1:end-4);
writecell(trial_trc, trial_fileName, 'Delimiter', '\t');
%
% Static
static_fileName = append(subjectName, trialID, '_static');
static = trial_trc(1:106,1:end);
static{1,4} = static_fileName;
static{3,3} = 100;
static{3,8} = 100;
writecell(static, static_fileName, 'Delimiter', '\t');

% Change files from .txt to .trc
file1 = append(trial_fileName, '.txt');
file2 = strrep(file1,'.txt','.trc');
copyfile(file1,file2);
%
file1 = append(static_fileName, '.txt');
file2 = strrep(file1,'.txt','.trc');
copyfile(file1,file2);


%% SET UP .MOT FILE

% Set up header info
trial_motName = append(subjectName, '_', trialID, '_GRF', '.mot');
dataCount = length(trial_data.time);

% Add specific header info
trial_mot{1,1} = trial_motName;
trial_mot{3,1} = append('nRows=',num2str(dataCount));

% Add frame, time, and marker coordinate data
for i = 1:dataCount
    trial_mot{i+7,1} = i*0.002-0.002; % Time
    trial_mot{i+7,2} = 0; % Left GRF x
    trial_mot{i+7,3} = rotated_data.GRF_LT(i); % Left GRF y
    trial_mot{i+7,4} = 0; % Left GRF z
    trial_mot{i+7,5} = rotated_data.CoP_LT_x(i); % Left CoP x
    trial_mot{i+7,6} = rotated_data.CoP_LT_y(i); % Left CoP y
    trial_mot{i+7,7} = rotated_data.CoP_LT_z(i); % Left CoP z
    %
    trial_mot{i+7,8} = 0; % Right GRF x
    trial_mot{i+7,9} = rotated_data.GRF_RT(i); % Right GRF y
    trial_mot{i+7,10} = 0; % Right GRF z
    trial_mot{i+7,11} = rotated_data.CoP_RT_x(i); % Right CoP x
    trial_mot{i+7,12} = rotated_data.CoP_RT_y(i); % Right CoP y
    trial_mot{i+7,13} = rotated_data.CoP_RT_z(i); % Right CoP z
    %
    trial_mot{i+7,14} = rotated_data.M_LT_x(i); % Left Moment x
    trial_mot{i+7,15} = 0; % Left Moment y
    trial_mot{i+7,16} = rotated_data.M_LT_z(i); % Left Moment z
    %
    trial_mot{i+7,17} = rotated_data.M_RT_x(i); % Right Moment x
    trial_mot{i+7,18} = 0; % Right Moment y
    trial_mot{i+7,19} = rotated_data.M_RT_z(i); % Right Moment z
end

% Remove <missing> in cells
mask = cellfun(@(x) any(isa(x,'missing')), trial_mot);
trial_mot(mask) = {[]};

%% WRITE GRF MOT FILES

trial_fileName = trial_motName(1:end-4);
writecell(trial_mot, trial_fileName, 'Delimiter', '\t');

% Change files from .txt to .trc
file1 = append(trial_fileName, '.txt');
file2 = strrep(file1,'.txt','.mot');
copyfile(file1,file2);