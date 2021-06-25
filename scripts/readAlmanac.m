%
% Author:    Batuhan UZUN
% Created:   23.06.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  ******** Week 115 almanac for PRN-01 ********
##  ID:                         01
##  Health:                     000
##  Eccentricity:               0.1102590561E-001
##  Time of Applicability(s):  319488.0000
##  Orbital Inclination(rad):   0.9847221889
##  Rate of Right Ascen(r/s):  -0.7737465154E-008
##  SQRT(A)  (m 1/2):           5153.642090
##  Right Ascen at Week(rad):   0.2345197091E+001
##  Argument of Perigee(rad):   0.863190763
##  Mean Anom(rad):             0.1332284785E+001
##  Af0(s):                     0.6494522095E-003
##  Af1(s/s):                  -0.1091393642E-010
##  week:                        115
##
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  - filename - enter the filename to be read.  If filename
%                      exists, the orbit will be calculated.
% 
% Output: - ALM - Output is a matrix with rows for each PRN and
%                       columns as follows:
% 
%   col  1:  id                   : satellite PRN of the SVN           
%   col  2:  health               : information about the state (health) of the 
%                                   entire GPS satellite constellation
%                                   (000=usable).
%   col  3:  eccentricity         : Amount of the orbit deviation from orbit.
%   col  4:  timeOfApplicability  : The number of seconds in the orbit when the 
%                                   almanac was generated.
%   col  5:  orbitalInclination   : The angle to which the SV orbit meets 
%                                   the equator.
%   col  6:  rateOfRightAscen     : Rate of change in the measurement of the 
%                                   angle of right ascension
%   col  7:  sqrtA                : The measurement from the center of the 
%                                   orbit to either the point of apogee or the 
%                                   point of perigee.
%   col  8:  rightAscenAtWeek     : Right Ascension at Time of Almanac (TOA)
%   col  9:  argumentOfPerigee    : An angular measurement measured from the 
%                                   ascending node to the point of perigee.
%   col 10:  meanAnomaly          : Angle (arc) traveled past the longitude 
%                                   of ascending node.
%   col 11:  af0                  : SV clock bias in seconds.
%   col 12:  af1                  : SV clock drift in seconds per seconds.
%   col 13:  week                 : GPS week (0000–1023), every 7 days since 
%                                   1999 August 22.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ALM, leapSecond] = readAlmanac( filename )

  almanacFile = fopen(filename);

  % is file opened ? 
  if almanacFile == -1
      errordlg(['The file ''' filename ''' does not exist.']);
      return;
  end

  i = 1;
  while feof(almanacFile) ~= 1

      % for title line
      currentLine = fgetl(almanacFile);

      %%% almanac contents %%%
      currentLine = fgetl(almanacFile);
      ALM.id(i) = str2num(currentLine(27:end));
      
      currentLine = fgetl(almanacFile);
      ALM.health(i) = str2num(currentLine(27:end));

      currentLine = fgetl(almanacFile);
      ALM.eccentricity(i) = str2double(currentLine(27:end));

      currentLine = fgetl(almanacFile);
      ALM.timeOfApplicability(i) = str2double(currentLine(27:end));

      currentLine = fgetl(almanacFile);
      ALM.orbitalInclination(i) = str2double(currentLine(27:end));

      currentLine = fgetl(almanacFile);
      ALM.rateOfRightAscen(i) = str2double(currentLine(27:end));

      currentLine = fgetl(almanacFile);
      ALM.sqrtA(i) = str2double(currentLine(27:end));

      currentLine = fgetl(almanacFile);
      ALM.rightAscenAtWeek(i) = str2double(currentLine(27:end));

      currentLine = fgetl(almanacFile);
      ALM.argumentOfPerigee(i) = str2double(currentLine(27:end));

      currentLine = fgetl(almanacFile);
      ALM.meanAnomaly(i) = str2double(currentLine(27:end));

      currentLine = fgetl(almanacFile);
      ALM.af0(i) = str2double(currentLine(27:end));

      currentLine = fgetl(almanacFile);
      ALM.af1(i) = str2double(currentLine(27:end));

      currentLine = fgetl(almanacFile);
      ALM.week(i) = str2num(currentLine(27:end));
      %%%------------------%%%
      
      % for empty line
      currentLine = fgetl(almanacFile);
      
      i = i+1;
  end
  
  % 2021 GPS-UTC Offset 
  leapSecond = 18; 
  
end