
# Naive X-Ray Simulator
Simulates x-ray images of one or more objects (STL files) created by an x-ray point source and a rectangular x-ray detector. The resulting simulation is visualized in a 3D plot and the resulting x-ray image written as a bitmap file.

X-ray attenuation is modeled according to exponential decay, i.e.:

I/I_0=exp(-Ax)

where I_0 and I are, respectively, the initial and attenuated x-ray intensities, x the path length through the material travelled and A the attenuation coefficient of the material. Values in the resulting image correspond to 1-(I/I_0), range of 0-1, inclusive.

## Getting Started and Usage
To install, unzip the supplied file and add the folder 'xraySimulator' to the path. All required functions/libraries are included and automatically added to the path when the simulation is run.

### Basic syntax:
```
image = xraySimulator('inputStlFilename.stl', attenuation, 'outputImageFilename.bmp');
```
Where attenuation is the x-ray attenuation (A in the above equation, units: cm^-1). Parameters of the simulation (object to source and object to detector distances, detector size and resolution, etc.) are set to the defaults defined in getDefaultParameters.m. Example:
```
img = xraySimulator('female_pelvis_fixed.stl',1, 'xray.bmp');
```
### Custom Parameters
 ```getDefaultParams()``` returns a structure containing the default values of various scene parameters. Parameters in the struct may be modified and the structure passed as an argument. Example:
```
myParams = getDefaultParams();
myParams.phantomToSourceDistance = 40; %Units - cm
img = xraySimulator('female_pelvis_fixed.stl', 1 , 'xray.bmp',myParams);
```
### Moving/Scaling Objects
Ojbects may be moved relative to the x-ray source/detector and/or scaled isotropically. Rotation (degrees), translation (cm) and scaling may be defined by a 1x7 vector:  ``` [Rx Ry Rz Tx Ty Tz S] ```.  Example (90 deg rotation, 5 cm translation in Z; scale x2):
```
img = xraySimulator('female_pelvis_fixed.stl', 1 , 'xray.bmp',myParams,[0 0 90 0 0 5 2]); % Custom parameters
img = xraySimulator('female_pelvis_fixed.stl', 1 , 'xray.bmp',[],[0 0 90 0 0 5 2]); % Default parameters
```
A scale value of 0 (non-existent object) is automatically corrected to 1 and a warning issued.

### Multiple Objects
A string array of N filenames and numeric array of N x-ray attenuations may be used to simulate multiple objects. Motion/scaling may still be specified as above using an Nx7 array. Example:
```
xraySimulator(["female_pelvis_fixed.stl"; "hipProsthesis.stl"], [0.1 0.5], 'xray.bmp',[],[0 0 -90 0 0 16.5 2; 0 0 90 -3 9 0 1]);
```
Output image of the above objects/positions with 512x512 detector size included as xray_512X512_hipAndProsthesis.bmp and xray_512X512_hipAndProsthesis.fig. Run time of 11m36s on MacBook Air (Early 2015), 1.6 GHz Intel Core i5 processor.

## Implementation Description
Intensity at each pixel is modeled by casting a line segment from the x-ray point source to the centre of each pixel. For each object (STL file), path length of the line segment through the object is calculated by finding the intersection points between the line segment and the object, sorting the points by proximity to the source and summing the distance between entrance and exit points (odd and even intersections, respectively). Path length is multiplied by the object's attenuation coefficient and returned as log-attenutation (i.e. A*x in previous equation).

The above process is repeated for all objects and the log-attenuation for each line segment summed. Relative intensity (I/I_0) is calculated and then transformed (1 - I/I_0) to generate pixel intensities.

X-ray/object intersections are found by testing for intersection of the line segment with all of the triangles in the object (mesh). This is time consuming due to the large number of triangles. For speed, a 3D rectangular bounding box is created around the object and intersection with the bounding box checked first (12 triangles in box vs thousands in typical object). Only line segments that intersect the bounding box are checked for intersection with the object. (For objects above, ~3-4x faster)

## Assumptions and Known Limitations
Higher order phenomenon (e.g. photon scattering, finite source size, detector noise) are not accounted for. Objects are assumed to be homogenous with respect to attenuation coefficient and attenuation coefficient is assumed to be constant for all components of the x-ray beam (e.g. mono-energetic beam or neglect of beam-hardening effects).

Simulations with multiple objects are not checked for object overlap/interference. Overlapping regions will appear to have an attenuation coefficient equal to the sum of the attenuation coefficients of all overlapping objects.

The currently implemented strategy relies on maintaining the surface representation of objects as described in the input STL files. While this approach is suitable for simulating well-delineated structures, it does not extend easily to use with volumetric data (e.g. CT scans) where structures are not explicitly outlined. While volumetric data could theoretically be simulated using this software by treating each voxel as a distinct structure with attenuation coefficient determined by the intensity at that voxel, it would be highly impractical. If voxelized data, rather than surfaces, was to be simulated, an alternate implementation of the sub-routine ```findXrayAttenuation()``` would need to be written.

Some additional features related to validating inputs and/or ensuring correct handling of edge cases have not yet been addressed. Where they have been identified, they are indicated in the source code comments with 'TODO' and a description.

## Acknowledgements

 - Intersections between objects (STL files) and line segments defining x-ray trajectories found using TriangleRayIntersections, written by Jaroslaw Tuszynski and published on Matlab Central File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection)
 - Collision detection between the x-ray source and the objects completed using in_polyhedron, written by Jaroslaw Tuszynski and published on Matlab Central File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/48041-in-polyhedron)
- Reading of STL files conducted using StlTools toolbox, written by Pau Micó and published on Matlab Central File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/51200-stltools)
- Display of progress bar during attenuation calculation completed using textprogressbar, written by Paul and published on Matlab Central File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/28067-text-progress-bar)
- STL file of hip prosthesis adapted from file from Edmundo Fuentes, published on GrabCad (https://grabcad.com/library/hip-replacement-prosthesis)
- STL file of pelvis adapted from file from Mike Barbee, published on Thingiverse (https://www.thingiverse.com/thing:1154961/#files)
