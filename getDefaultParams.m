function [ dp ] = getDefaultParams()
%GETDEFAULTSCENEPARAMS Summary of this function goes here
%   Detailed explanation goes here

%TODO Generates default parameters for x-ray simulation that are
%independent of the phantom being imaged (e.g. detector size, source location,
%etc.). Also creates template to be modified by user for custom parameters.


%TODO If time, convert params from structure to object to prevent
%mismatched data types, misnamed vars, etc. If done, also add check of
%object type to checkInputs


%Geometric/Simulation Parameters
dp.phantomToSourceDistance = 100; %cm
dp.phantomToDetectorDistance = 20; %cm
dp.detectorWidthDist= 40; %cm
dp.detectorWidthPixels= 32; 
dp.detectorHeightDist= 40; %cm
dp.detectorHeightPixels = 32;

%Other
dp.showOutput = true;
dp.phantomTransparency = 0.5;
dp.writeImageToDisk = true;

end

