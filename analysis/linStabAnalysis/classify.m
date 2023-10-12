function binNumber = classify(z, eigVec)
%   classify(numberOfBins, eigVec)
%       Determines whether a given eigenvector is a surface mode,
%       halocline mode, deep mode, etc. Specifically, defines normalized
%       vectors with non-zero components at associated depths and takes
%       inner product between that and eigVec. eigVec then gets
%       assigned an int (binNumber) marking how its been classified.
%
%
%       Inputs
%       levBounds   --  depths marking the boundary between sections of
%                       the water column. Two values => three bins, etc.
%       eigVec      --  the eigenvector to be classified; calculated in
%                       qggrz.m
%       z           --  cell centers of vertical grid
%
%       Output
%       binNumber   --  an int ranging from 1 to length(levBounds)+1 that
%                   --  gives the classification for eigVec. 1 is surface
%       
%       Params
%       lev1        --  depth (m) between surface and halocline section
%       lev2        --  depth (m) between halocline and deep section

% params
lev1    = 75.1;
lev2    = 250.1;
numLevs = length(z);

% initialize vecs
vec1 = zeros(numLevs);
vec2 = zeros(numLevs);
vec3 = zeros(numLevs);

for index = numLevs
    if lev1 > z(index)
        vec1(index) = 1;
    end
    if (lev1 < z(index)) && (z(index) < lev2)
        vec2(index) = 1;
    end
    if (lev2 < z(index))
        vec3(index) = 1;
    end
end

% normalize vecs
vec1 = vec1 / dot(vec1,vec1);
vec2 = vec2 / dot(vec2,vec2);
vec3 = vec3 / dot(vec3,vec3);

% calculate inner producs
val1 = dot(vec1,eigVec);
val2 = dot(vec2,eigVec);
val3 = dot(vec3,eigVec);

% classify eigVec
if max([val1 val2 val3]) == val1
    binNumber = 1;
elseif max([val1 val2 val3]) == val2
    binNumber = 2;
elseif max([val1 val2 val3]) == val3
    binNumber = 3;
else
    binNumber = 0;
end

% return value
