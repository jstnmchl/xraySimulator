function [ image ] = xraySimulator(stlFilename, attenCoeffs, outputImgFilename, params, varargin ) 
%XRAYSIMULATOR Simulates x-ray images of one or more objects (STL files) 
%   created by an x-ray point source and a rectangular x-ray detector. The 
%   resulting simulation is visualized in a 3D plot and the simulated x-ray 
%   image written as a bitmap file. 
%   
%   For detailed description and usage, see xraySimulator_README.md

tic

if (nargin<4 || isempty(params) ) 
    params = getDefaultParams();
end

if nargin > 4
    rotTransScale = varargin{1};
else
    rotTransScale =[];
end

if nargin > 5
    error('Too many input variables')
end

addLibsToPath();
checkInputs();

phantom = generateVirtualPhantom(stlFilename, attenCoeffs, rotTransScale);
scene = generateSimulationScene(params);
image = generateSimulatedImage(phantom, scene);

if params.showOutput
    plotScene(scene,phantom,image,params)
end

if params.writeImageToDisk
    writeImgToDisk(image,outputImgFilename);
end

runTime = toc;
disp(['Total run time: ' num2str(runTime) ' seconds'])

    function [] = checkInputs()
    %%% Check validity of inputs
        
    %Ensure stlFilename is string array (if only one filename, convert
    %char to string). 
    if not(isstring(stlFilename))
        if ischar(stlFilename) && (size(stlFilename,1) == 1)
            stlFilename = string(stlFilename);
        else
            error('Invalid variable type of stlFilename. Variable should be string array or 1D character array');
        end
    end

    %TODO: does stlfilename end in .stl?

    %TODO: does imgFilename end in .bmp?

    %TODO: Same number of stls as attenuation coeffs?

    %TODO: check distances in params are non-negative

    %TODO: check pixels in params (width/height) >0

    %TODO: Check all variables in params defined, no other params defined (to
    %protect against misspelled variable names)
    
    %If present, check rotTransScale for correct size, non-negative
    %scale. If scale = 0, change to 1 and display warning.
    if ~isempty(rotTransScale)
        for i=1:size(rotTransScale,1)
            if rotTransScale(i,7) < 0
                error(['Scale value must be non-negative. Scale of ' char(stlFilename(i)) ' was ' num2str(rotTransScale(i,7)) ])
            elseif rotTransScale(i,7) == 0
                rotTransScale(i,7) = 1;
                warning(['Scale value for ' char(stlFilename(i)) ' was 0 (implies non-existent object). Automatically corrected to 1.'])
            end
        end
    end

      

    end
end

function [] = addLibsToPath()
%Adds folders of third-party functions to path. Temporarily changes CD to location of 
%xraySimulator to ensure folders can be found

scriptFolder = fileparts(mfilename('fullpath'));

oldFolder = cd(scriptFolder);

addpath(genpath('3rdParty')) 

cd(oldFolder);

end


function [p] = generateVirtualPhantom(stlFilename, attenCoeffs, rotTransScale) 
%Reads stl(s) from file and save to struct of faces, vertices and filenam.
%Associate X-ray atenuation coefficient(s) with each phantom. 
numStls = size(stlFilename,1);

p = cell(numStls,1); %phantom

for i=1:numStls
    [p{i}.vertices,p{i}.faces,~,~] = stlRead(stlFilename(i));
    p{i}.name = stlFilename(i);
    p{i}.attenCoeff = attenCoeffs(i);
    
    if ~isempty(rotTransScale)
        p{i}.vertices = transformVertices(p{i}.vertices,rotTransScale(i,:));
    end
end

end

function v = transformVertices(v, rts)
%Transforms vertices of object according to rts (rotation - deg, translation - cm, 
%scale). Rotation and scaling applied about object centre (i.e. mean of vertices).
%rts = [Rx Ry Rz Tx Ty Tz S];

%centre
centre = mean(v);
v = bsxfun(@minus,v,centre);

%rotate
dcm = angle2dcm(deg2rad(rts(3)),deg2rad(rts(2)),deg2rad(rts(1)), 'ZYX');
v = v*dcm;

%scale
scale = eye(3).*rts(7);
v = v*scale;

%un-centre
v = bsxfun(@plus,v,centre);

%translate
v = bsxfun(@plus,v,rts(4:6));
end


function [scene] = generateSimulationScene(params)
%Generates geometry of scene (source, detector, rays) in
%common coordinate system with phantom.

%Origin of scene corresponds to origin of STL file. Source lies on X axis
%(x-coordinate < 0). Detector is perpendicular to x axis (x-coordinate > 0),
%centred on x-axis.

scene.sourceCoords = [-1*params.phantomToSourceDistance 0 0];

scene.pixelCoords = generateDetectorPixelCoords(params);

scene.lineSegs = generateSourceToDetectorLineSegs(scene.sourceCoords,scene.pixelCoords);
end

function [pixels] = generateDetectorPixelCoords(params)
%Generate 3xNxM array defining centres of NxM pixels of detector

%numPixels = params.detectorWidthPixels*params.detectorHeightPixels;
%pixels = zeros(numPixels,3);

pixels = zeros(3, params.detectorWidthPixels,params.detectorHeightPixels);

pixelWidth = params.detectorWidthDist/params.detectorWidthPixels;
pixelHeight = params.detectorHeightDist/params.detectorHeightPixels;

xVal = params.phantomToDetectorDistance;
yVals = -1*(0.5*params.detectorWidthDist-0.5*pixelWidth):pixelWidth:(0.5*params.detectorWidthDist-0.5*pixelWidth);
zVals = -1*(0.5*params.detectorHeightDist-0.5*pixelHeight):pixelHeight:(0.5*params.detectorHeightDist-0.5*pixelHeight);

[Y,Z] = meshgrid(yVals, zVals);
pixels(1,:,:) = repmat(xVal,params.detectorHeightPixels,params.detectorWidthPixels);
pixels(2,:,:) = Y;
pixels(3,:,:) = Z;
end

function [lineSegs] = generateSourceToDetectorLineSegs(source, pixels)
%Returns 3xNxM array, defining line segments from the source to 
%each of NxM pixels as vectors

lineSegs=zeros(size(pixels));

[~, height, width] = size(pixels);

for h=1:height
    for w=1:width
        lineSegs(:,h,w) = pixels(:,h,w)' - source;
    end
end

end


function [img] = generateSimulatedImage(phantom, scene)
%Calculates simulated Xray image from scene, phantom geometry

checkForCollisions(phantom,scene);

numObjects = size(phantom,1);

height = size(scene.lineSegs,2);
width = size(scene.lineSegs,3);
attenuation = zeros([height width numObjects]);

for i=1:numObjects
    attenuation(:,:,i) = findXrayAttenuation(phantom{i}, scene);
end

%Find relative intensity (i.e. I/I_0) at each pixel
rIntensity = exp(-1*sum(attenuation,3));

%Invert rIntensity values so dense objects (e.g. bone) appear bright
img = -1*rIntensity + 1;

%Correct for direction conventions
img = flipud(fliplr(img));
end

function [] = checkForCollisions(phantom,scene)
%Checks for collisions between objects and source. (Collisions between
%objects not checked for, treated as summative overlap)

%TODO Check for collisions between objects and detector

numObjects = size(phantom,1);
isSourceInside = zeros(numObjects,1);
for i=1:numObjects
    isSourceInside(i) = in_polyhedron(phantom{i}.faces, phantom{i}.vertices,scene.sourceCoords);
end

if any(isSourceInside)
    names = [];
    for i=1:numObjects
        if isSourceInside(i)
            names = [names ' ' char(phantom{i}.name)];
        end
    end
    error([num2str(sum(isSourceInside)) ' STL files collide with the x-ray source. Colliding files are:' names ])
end

end

function [atten] = findXrayAttenuation(object, scene)
%Returns NxM array of attenuation (i.e. x in I=I_0*e^-x )for a single 
%closed surface where NxM is the number of pixels (width x height) on the detector 

boundingBox = generateBoundingBox(object);
%boundingBox=[];

[~, height, width] = size(scene.pixelCoords);

atten=zeros(height,width);

%Create progress bar for finding attenuation (can be lengthy)
numProjections = numel(atten);
progress=0;

textprogressbar(['Projecting X-Rays for ' char(object.name) ': ']);
%If progress bar interrupted, previous bar will not be closed and result in 
%error on subsequent run. Calling textprogressbar a second time fixes this
%problem.
try 
    textprogressbar(0);
catch
    textprogressbar('');
end

for h=1:height
    for w=1:width
        atten(h,w) = estLineIntegral(scene.sourceCoords, scene.lineSegs(:,h,w), object, boundingBox);
        progress = progress+1;
        textprogressbar(progress/numProjections *100);
    end
end

textprogressbar(' Done.');

end

function [box] = generateBoundingBox(obj)
%Returns vertices of rectangular bounding box around object in format
%required by TriangleRayIntersections
v = obj.vertices;
boxRange =[min(v(:,1)) max(v(:,1)); min(v(:,2)) max(v(:,2)); min(v(:,3)) max(v(:,3))];
[X, Y, Z] = meshgrid(boxRange(1,:),boxRange(2,:), boxRange(3,:));

boxVertices = [reshape(X,[],1,1) reshape(Y,[],1,1) reshape(Z,[],1,1)];

boxFaces = boundary(boxVertices(:,1),boxVertices(:,2),boxVertices(:,3));

box(:,:,1) = boxVertices(boxFaces(:,1),:);
box(:,:,2) = boxVertices(boxFaces(:,2),:);
box(:,:,3) = boxVertices(boxFaces(:,3),:);

end

function [integral] = estLineIntegral(src,lineSeg, object, bBox)
%Returns estimate of line integral from source to centre of one
%detector pixel through object by finding path length through object,
%multiplying by object's x-ray attenuation.

%Check intersection with bounding box first for speed
[intersect_box, ~, ~, ~,~] = TriangleRayIntersection(src, lineSeg, bBox(:,:,1), bBox(:,:,2), bBox(:,:,3),'lineType','segment');

if sum(intersect_box) == 0
    integral = 0;
    return
end

vert1 = object.vertices(object.faces(:,1),:);
vert2 = object.vertices(object.faces(:,2),:);
vert3 = object.vertices(object.faces(:,3),:);

[intersect, T, ~, ~,xcoor] = TriangleRayIntersection(src, lineSeg, vert1, vert2, vert3,'lineType','segment');

numIntersections = sum(intersect);

pathLength = 0;

if numIntersections == 0
    integral = 0;
    return
elseif mod(numIntersections,2)
    %TODO Handle edge case (pun intended) where ray intersects with edge of 
    %outermost triangle, leading to only one intersection (i.e. entrance 
    %and exit are same point). Should return integral = 0
    error(['Odd number of intersections. X-ray appears to enter object but not exit or vice versa. Object: ' char(object.name)] )
end

intersectCoords = xcoor(intersect,:);

%Sort intersections by proximity to source
intersectT = T(intersect);
intersectCoords = [intersectCoords intersectT];
intersectCoords = sortrows(intersectCoords,4);
intersectCoords = intersectCoords(:,1:3);


%Assuming odd intersections are x-ray entering object & even leaving, sum
%path length through object
for i=2:2:numIntersections
    pathLength = pathLength + norm(intersectCoords(i,:) - intersectCoords(i-1,:));
end

integral = pathLength * object.attenCoeff;
end


function [] = plotScene(scene,phantom, img, params)

figure()
%Plot source
plot3(scene.sourceCoords(1),scene.sourceCoords(2),scene.sourceCoords(3),'r*');

grid on; hold on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

%Plot detector/image
W = params.detectorWidthDist/2;
H = params.detectorHeightDist/2;

X = repmat(params.phantomToDetectorDistance,2,2);
Y = [W -W; W -W];
Z = [H H; -H -H];

surface(X,Y,Z,repmat(img,1,1,3), 'FaceColor','texturemap','EdgeColor','none');

numObjects = size(phantom,1);
C = {'b','r','y',[.5 .6 .7],[.8 .2 .6], 'k','g'}; % Cell array of colors.
for i=1:numObjects
    cInd = mod(i,numel(C)); %Cycle through colors
    obj.vertices = phantom{i}.vertices;
    obj.faces = phantom{i}.faces;
    patch(obj,'FaceColor', C{cInd}, ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15, ...
         'FaceAlpha', params.phantomTransparency);
end


end


function [] = writeImgToDisk(img,outputImgFilename)

%TODO handle scenario where filename already exists

imwrite(img, outputImgFilename);

end