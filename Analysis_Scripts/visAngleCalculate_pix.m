% calculate visual angle based on object width (cm), distance from screen
% (cm) and the screen dimensions (diagonal, inches)

function thisVisAngle = visAngleCalculate_pix(objectWidth_pix,screenDist_cm)

% assuming 1280 * 1024 resolution on a 19 inch monitor, each pixel = .2944
% mm = .02944 cm

% ACTUAL screen measurement for EEG chamber is 18 inch monitor, with width
% = 34.8 (actual viewing width) ...at 1280 x 1024 res this = 34.8/1024 = 

%pixelWidth = .02944; % cm
%pixelWidth = .0482; %???
pixelWidth = .03352; % ORIGINAL


objectWidth_cm = objectWidth_pix * pixelWidth; % convert object width to cm

thisVisAngle = 2*atan(objectWidth_cm/(2*screenDist_cm))*180/pi; % calculate visual angle

end