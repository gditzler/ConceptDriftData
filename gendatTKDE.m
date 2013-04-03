%% About this file
%  This file is meant to sevre as a demo on how to call the functions
%  which generate the data presented in [1]. However, in [1], the data
%  are drawn from an unbalanced source. Therefore, one of the classes
%  must be undersampled to created an environment that contains concept
%  drift and class imbalance. Class '2' was undersampled for every 
%  dataset tested in [1].
%  
%  References
%   1. Ditzler, G. & Polikar, R
%      "Incremental Learning of Concept Drift from Streaming Imbalanced Data"
%      Submitted to IEEE Transactions on Knowledge and Data Engineering, 2012
%      accepted. 
%     
%
%  gendatTKDE.m
%  Copyright (C) 2011  Gregory Ditzler and Robi Polikar
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%% Gaussian Data
clc;clear;close all;
disp('Generating drift Gaussian data. Press')
disp('any key to continue.')
disp(' ');pause

T = 100;
N = 500;
[xTrain,cTrain,xTest,cTest] = ConceptDriftData('gaussian',T,N);
figure;
for t = 1:T
  clf;hold on;
  plot(xTrain{t}(1,cTrain{t}==1),xTrain{t}(2,cTrain{t}==1),'r.');
  plot(xTrain{t}(1,cTrain{t}==2),xTrain{t}(2,cTrain{t}==2),'g.');
  axis([-2,10,-2,10])
  pause(0.1);
end
%% SEA Data
clc;clear;close all;
disp('Generating drift SEA data. Press')
disp('any key to continue.')
disp(' ');pause

T = 200;
N = 1000;
[xTrain,cTrain,xTest,cTest] = ConceptDriftData('sea',T,N);
figure;
for t = 1:T
  clf;hold on;
  plot(xTrain{t}(1,cTrain{t}==1),xTrain{t}(2,cTrain{t}==1),'r.');
  plot(xTrain{t}(1,cTrain{t}==2),xTrain{t}(2,cTrain{t}==2),'g.');
  axis([0,10,0,10])
  pause(0.1);
end
%% CB Data
clc;clear;close all;
disp('Generating drift Checkerboard data. Press')
disp('any key to continue.')
disp(' ');pause

T = 300;
N = 2000;
[xTrain,cTrain,xTest,cTest] = ConceptDriftData('checkerboard',T,N);
figure;
for t = 1:T
  clf;hold on;
  plot(xTrain{t}(1,cTrain{t}==1),xTrain{t}(2,cTrain{t}==1),'r.');
  plot(xTrain{t}(1,cTrain{t}==2),xTrain{t}(2,cTrain{t}==2),'g.');
  axis([0,1,0,1])
  pause(0.1);
end
%% SP Data
clc;clear;close all;
disp('Generating drift Checkerboard data. Press')
disp('any key to continue.')
disp(' ');pause

T = 300;
N = 500;
[xTrain,cTrain,xTest,cTest] = ConceptDriftData('spiral',T,N);
figure;
for t = 1:T
  clf;hold on;
  plot(xTrain{t}(1,cTrain{t}==1),xTrain{t}(2,cTrain{t}==1),'r.');
  plot(xTrain{t}(1,cTrain{t}==2),xTrain{t}(2,cTrain{t}==2),'g.');
  axis([-2*pi,2*pi,-2*pi,2*pi])
  pause(0.1);
end