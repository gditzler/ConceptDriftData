% CCONCEPTDRIFTDATA: Generate Data Streams with Concept Drift
%   [xTrain,cTrain,xTest,cTest] = conceptDriftData(type,T,N)
% 
%   This m-file generates data containing concept drift and
%   classes can be under-sampled to create an imbalanced 
%   learning scenario. The rotating checkerboard function was
%   written by Ludmila Kuncheva.
%   
%   GENDATCB link
%   http://pages.bangor.ac.uk/~mas00a/data_public/artificial_data/gendatcb.m
%   
%   INPUT
%     TYPE: Dataset desired
%       'noaa': Weather dataset collected from the NOAA's website. 
%         The data contains 8 features and approximately 50 years
%         worth of data. 
%       'elec2': VIC/NSW Electricity pricing dataset used by 
%         several researchers involved with concept drift. Two 
%         features have been removed as described in previous
%         papers.
%       'sea': Street & Kim's shifting hyperplane problem with
%         a small modification (4,7,4,7) rather than (7,9,8,9.5)
%       'checkerboard': Roatating checkerboard problem with a 
%         linear drift. The CB rotates over a range of 0 - 2pi
%         with 0-pi being a recurring environment at pi-2pi
%       'spiral': Generate a rotating spiral drift scenario that
%         contains two-class where each class contains two spirals.
%         This drift scenario contains reoccurring environments.
%       'gaussian': Guassian data is generated from 4 components
%         where three are majority class components and one is a
%         minority class components. Drift is simulated by 
%         varying the covariances and means.
%     T: Number of time stamps (ingored in 'elec2' and 'noaa')
%     N: Number of instances to generate
%   OUTPUT
%     xTrain,xTest: cell array containing the feature vectors in
%       column format for training/testing. Each entry to the 
%       cell array is an index to a time stamp
%     cTrain,cTest: cell array containing the instance labels in
%       corresponding to xTrain/xTest. Each entry to the 
%       cell array is an index to a time stamp
%

%
%     ConceptDriftData.m
%     Copyright (C) 2011  Gregory Ditzler, Ryan Elwell, and Robi Polikar
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
function [xTrain,cTrain,xTest,cTest] = ConceptDriftData(type,T,N)
  % switch on the type of experiment
  switch type
    case 'noaa'
      % T is disregarded; we will use all the data, so the 
      % size of the N will determine T.
      alldata=load('noaa_data.dat')';
      allclass=load('noaa_label.dat')';
      T = floor(floor(size(allclass,2)-N-1)/N)-1;
      for n=1:T
        starti = ((n-1)*N)+1;
        stopi = (starti+N);
        trinds = starti:stopi;
        tsinds = stopi+1:stopi+N;
        xTrain{n} = alldata(:,trinds);
        cTrain{n} = allclass(:,trinds);
        xTest{n} = alldata(:,tsinds);
        cTest{n} = allclass(:,tsinds);
      end    
    case 'elec2'
      alldata=load('elec2_data.dat')';
      allclass=load('elec2_label.dat')';
      alldata(3,:)=[]; % remove price features
      alldata(4,:)=[];
      T = floor(floor(size(allclass,2)-N-1)/N)-1;
      for n=1:T
        starti = ((n-1)*N)+1;
        stopi = (starti+N);
        trinds = starti:stopi;
        tsinds = stopi+1:stopi+N;
        xTrain{n} = alldata(:,trinds);
        cTrain{n} = allclass(:,trinds);
        xTest{n} = alldata(:,tsinds);
        cTest{n} = allclass(:,tsinds);
      end
    case 'sea'
      theta = [4*ones(1,50),7*ones(1,50)];
      theta = [theta,theta];
      noise = 0.05;
     
      for t = 1:T
        [xTrain{t},cTrain{t}] = SEA(theta(t),N,noise);
        [xTest{t},cTest{t}] = SEA(theta(t),N,noise);
      end
    case 'checkerboard'
      a = linspace(0,2*pi,T);
      side = 0.5;
      
      for t = 1:T
        [xTrain{t},cTrain{t}] = CBDAT(side,a(t),N);
        [xTest{t},cTest{t}] = CBDAT(side,a(t),N);
      end
    case 'spiral'
      range = 2*pi;
      spirals = 4;
      t0=linspace(0,range,N);
     
      for t = 1:T
        phi=2*pi*(t-1)/T;
        TS=[cos(phi) sin(phi);
           -sin(phi) cos(phi)];
       
        [X1,Y1] = SPDAT(t0,TS,spirals);
        [X2,Y2] = SPDAT(t0,TS,spirals);

        xTrain{t} = X1;
        cTrain{t} = Y1;
        xTest{t} = X2;
        cTest{t} = Y2;
      end
    case 'gaussian'
      t1 = 1/3;
      t2 = 2/3;
      m = 1;
      
      for t = 0:1/T:1-1/T
        if (t>=0)&&(t<=t1),
          mu1 = [2 5];
          mu2 = [8 5];
          mu3 = [5 2];
          mu4 = [5 8];
          sigma1 = [1 0;0 (1+6*t)];
          sigma2 = [1 0;0 1];
          sigma3 = [(3-6*t) 0;0 1];
          sigma4 = [(3-6*t) 0;0 1];

          x1 = mvnrnd(mu1,sigma1,2*N)';
          x2 = mvnrnd(mu2,sigma2,2*N)';
          x3 = mvnrnd(mu3,sigma3,2*N)';
          x4 = mvnrnd(mu4,sigma4,2*N)';
        elseif (t>t1)&&(t<=t2),
          mu1 = [2 5];
          mu2 = [(8-9*(t-t1)) 5];
          mu3 = [(5+9*(t-t1)) 2];
          mu4 = [(5+9*(t-t1)) 8];
          sigma1 = [1 0;0 3];
          sigma2 = [1 0;0 1];
          sigma3 = [1 0;0 1];
          sigma4 = [1 0;0 1];

          x1 = mvnrnd(mu1,sigma1,2*N)';
          x2 = mvnrnd(mu2,sigma2,2*N)';
          x3 = mvnrnd(mu3,sigma3,2*N)';
          x4 = mvnrnd(mu4,sigma4,2*N)';
        else
          mu1 = [(2+6*(t-t2)) (5-9*(t-t2))];
          mu2 = [(5-3*(t-t2)) (5+9*(t-t2))];
          mu3 = [8 2];
          mu4 = [8 8];
          sigma1 = [1 0;0 (3-6*(t-t2))];
          sigma2 = [1 0;0 1];
          sigma3 = [1 0;0 1];
          sigma4 = [1 0;0 1];

          x1 = mvnrnd(mu1,sigma1,2*N)';
          x2 = mvnrnd(mu2,sigma2,2*N)';
          x3 = mvnrnd(mu3,sigma3,2*N)';
          x4 = mvnrnd(mu4,sigma4,2*N)';
        end
        xTrain{m} = [x1(:,1:N),x2(:,1:N),...
          x3(:,1:N),x4(:,1:N)];
        cTrain{m} = [1*ones(1,N),2*ones(1,N),...
          1*ones(1,N),1*ones(1,N)];
        
        xTest{m} = [x1(:,N+1:end),x2(:,N+1:end),...
          x3(:,N+1:end),x4(:,N+1:end)];
        cTest{m} = [1*ones(1,numel(N+1:2*N)),2*ones(1,numel(N+1:2*N)),...
          1*ones(1,numel(N+1:2*N)),1*ones(1,numel(N+1:2*N))];
        m = m+1;
      end
    otherwise
      error('ERROR::conceptDriftData.m: Unknow Dataset Selected');
  end
end
%% extra functions
function [X,Y] = SEA(theta,N,noise)
  Nnoise = floor(N*noise);
  %%%% generate data and labels associated with each instance
  x = 10*rand(3,N);
  X = x; % save data in X
  Y = zeros(1,length(X));
  x(3,:) = []; % remove feature w/ no information
  Y(sum(x)>=theta) = 1;
  Y(Y~=1) = 2;
  %%%% add noise into the dataset
  r = randperm(numel(Y));
  r(Nnoise+1:numel(Y)) = [];
  Y(r) = Y(r)-1; % change 2->1 and 1->0
  Y(Y==0) = 2;   % class '0' is actually class '2' with noise
end

function [X,Y] = CBDAT(a,alpha,N)
  % L. I. Kuncheva, "Combining Pattern Classifiers", Wiley & Sons, 2004
  [X,Y] = gendatcb(N,a,alpha);
  X = X';
  Y = Y';
end

function [X,Y] = SPDAT(t0,TS,spirals)
  angle = 2*pi/spirals;
  X = [];Y=X;
  b = 0.5;% control the noise in the spirals
  for i = 1:spirals
    [x,y] = parametric(t0);
    sc = cos((i-1)*angle);  
    ss = sin((i-1)*angle);
    % rotation transformation matrix
    T=[ sc ss; 
        ss  -sc];
    %combine all the components 
    A=[x; y]';  

    Z = A*T; 
    Z = Z*TS;

    X=[X,Z'];
    if mod(i,2)==0,
      Y=[Y,2*ones(1,numel(t0))];
    else
      Y=[Y,1*ones(1,numel(t0))];
    end
  end
  X = X+b*randn(size(X,1),size(X,2));
end

function [x y] = parametric(theta)
  x = theta.*cos(theta);
  y = theta.*sin(theta);
end

function [d,labd] = gendatcb(N,a,alpha)
% N data points, uniform distribution,
% checkerboard with side a, rotated at alpha
d = rand(N,2);
d_transformed = [d(:,1)*cos(alpha)-d(:,2)*sin(alpha), ...
    d(:,1)*sin(alpha)+d(:,2)*cos(alpha)];
s = ceil(d_transformed(:,1)/a)+floor(d_transformed(:,2)/a);
labd = 2-mod(s,2);
end