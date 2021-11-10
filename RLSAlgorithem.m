clc
clear all
close all
%% RECURSIVE LEAST SQUARES ALGORITHEM.
%   This program was developed as a course project for Fuzzy systems course by
%    Habibollah Naeimi.
disp(' RECURSIVE LEAST SQUARES ALGORITHEM.');
disp('  This program was developed as a course project for Fuzzy systems course')
disp('   by Habibollah Naeimi.')
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 1st Part: Parameter Setting.
%%  Parameters Initiating.
disp(' Parameters Initiating...');

Sigma = 100;                       
%Sigma = input('Please Enter Sigma(The Larg Constant!):\n');

                                       % Sampling paramers
DataPairsNum = 250;                    % Number of Data Pairs.
%DataPairsNum = input('Please Enter Number of Data Pairs:\n');
SamplesNum = 500;                      % Number of Samples.
%SamplesNumb = input('Please Enter Number of Samples:\n');
MFN = 6;                               % Number of Membership Function.
%MFN = input('Please Enter Number of Membership Function:\n');
MFType = 3;                            % Type of Membership Function.
%MFType = input('Please Enter Type of Membership Function:\nTri=1  Trap=2  Gauss=3 MFType:');

InpNum = 4;                            % Inputs Number.
LowBnd = 0.2;                          % Lower Bound.
UpBnd = 1.4;                           % Upper Bound.

%%  Fixing Parameters dimentions.

e1 = numel(MFN);
e2 = numel(LowBnd);
e3 = numel(UpBnd);
e4 = numel(MFType);

if e1~=InpNum+1
    if e1==1                           % Fixing MF Number.
        MFN = repmat(MFN,1,InpNum+1);
    else
        disp(' Invalid Membership Function Number!');
        disp(' Please Start Again.');
        MFN = nan;
    end
end

if e2~=InpNum+1
    if e2==1                           % Fixing Lower Bound.
        LowBnd = repmat(LowBnd,1,InpNum+1);
    else
        disp(' Invalid Lower Boundary!');
        disp(' Please Start Again.');
        LowBnd = nan;
    end
end

if e3~=InpNum+1
    if e3==1                           % Fixing Upper Bound.
        UpBnd = repmat(UpBnd,1,InpNum+1);
    else
        disp(' Invalid Upper Boundary!');
        disp(' Please Start Again.');
        UpBnd = nan;
    end
end

if e4~=InpNum+1
    if e4==1                          % Fixing MF Type.
        MFType = repmat(MFType,1,InpNum+1);
    else
        disp(' Invalid Membership Function Type!');
        disp(' Please Start Again.');
        MFType = nan;
    end
end

disp(' Part 1: DONE!');
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 2nd Part: Sampling.
%%  Calculating Samples.

SAMPLES = zeros(SamplesNum,InpNum+1);
Samples1 = 0.2:0.01:0.51;

for i=33:SamplesNum+33+InpNum
    Samples1(i) = 0.2*Samples1(i-31)/(1+(Samples1(i-31)^10))+0.9*Samples1(i-1);
end

Samples1 = Samples1(33:end);
    
for i=1:SamplesNum
    SAMPLES(i,:) = Samples1(i:i+InpNum);
end

Pairs = SAMPLES(1:DataPairsNum,:);
    
disp(' Time Series Sampling is Reasdy!');
disp(' ');
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 3rd Part: Rule Base.
%%  Rule base.

Rules = RuleBase(MFN(1:end-1),InpNum);

RulesNumber = numel(Rules(:,1));
disp(' Number of Rules is:');
disp(RulesNumber);
disp(' Rule Base is Reasdy!');
disp(' ');
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 4th Part: Main Computations.
%%  Output Center Generation.
% Sth like online initial parameters choosing.
Theta = LowBnd(end):(UpBnd(end)-LowBnd(end))/(size(Rules,1)-1):UpBnd(end);
%Theta=rand(1,size(Rules,1))*(UpperBound(end)-LowerBound(end))+LowerBound(end);
Theta = Theta';                            % Fixing Dimentions.

%%  RLS Computations.
P = Sigma*eye(size(Rules,1));

for p=1:size(Pairs,1)
    
    b_x = CalculatingB(Pairs(p,1:end-1),Rules,MFN,MFType,UpBnd,LowBnd);
    K_p = P*b_x*(1/(b_x'*P*b_x+1));
    Theta = Theta+K_p*(Pairs(p,end)-b_x'*Theta);
    P = P-P*b_x*(1/(b_x'*P*b_x+1))*b_x'*P;
    
end

disp(' Primary Computations are done.');
disp('*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*-.-*------------------------------------------------------')
disp(' ');
%% 5th Part: Results.
%%  Results Calculations.
f = zeros(1,SamplesNum);
f(1:2) = SAMPLES(1:2,end)';

for i=2:size(SAMPLES,1)
    
    b = CalculatingB(SAMPLES(i,1:end-1),Rules,MFN,MFType,UpBnd,LowBnd);
    f(i) = b'*Theta;
    y_aprx(i) = f(i);
               
    
end

%%  Plotting Results.
figure;
plot(SAMPLES(:,end));
hold on
plot(y_aprx,'r');
legend('Real Vaue','Predicted');

%%  Errors.
disp('Mean Square Error:')
MSE = mse(SAMPLES(:,end)-y_aprx')
disp('Mean Absolute Error:')
MAE = mae(SAMPLES(:,end)-y_aprx')
