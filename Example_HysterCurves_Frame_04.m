% Example_HysterCurves_Frame_04
%----------------------------------------------------------------
% PURPOSE 
%    To compute the hysteresis curves for a
%    reinforced concrete plane frame
%
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-02-23
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------
clc
clear all

nnodes=10;
nbars=9;

%% Materials
% f'c of each element
fpc=[300;
    300;
    250;
    250;
    250;
    300;
    250;
    250;
    250];

% Elasticity modulus of each element in function of f'c
E=zeros(nbars,1);
for i=1:nbars
    E(i)=14000*sqrt(fpc(i));
end

%% Geometry/Topology
dimensions=[40 40;
            30 60;
            40 40;
            30 60;
            40 40;
            30 60;
            40 40;
            30 60;
            40 40];
        
% cross-section area of each element
A=zeros(nbars,1);
for i=1:nbars
    A(i)=dimensions(i,1)*dimensions(i,2);
end

I=zeros(nbars,1);
for i=1:nbars
    I(i)=1/12*dimensions(i,1)*dimensions(i,2)^3;
end

% coordinates of each node for each bar
coordxy=[0 0;
         0 300;
         500 0;
         500 300;
         800 0;
         800 300;
         1400 0;
         1400 300;
         1600 0;
         1600 300]; 

% final and initial nodes for each element
ni=[1;2;3;4;5;6;7;8; 9];
nf=[2;4;4;6;6;8;8;10;10];

% prescribed boudnary conditions [dof, displacement]
bc=[1 0;
    2 0;
    3 0;
    7 0;
    8 0;
    9 0;
    13 0;
    14 0;
    15 0;
    19 0;
    20 0;
    21 0;
    25 0;
    26 0;
    27 0];

supports=[1 "Fixed" "Fixed";
           2 "Fixed" "Fixed";
           3 "Fixed" "Fixed";
           4 "Fixed" "Fixed";
           5 "Fixed" "Fixed";
           6 "Fixed" "Fixed";
           7 "Fixed" "Fixed";
           8 "Fixed" "Fixed";
           9 "Fixed" "Fixed"];
   
type_elem=[1 "Col";
           2 "Beam";
           3 "Col";
           4 "Beam";
           5 "Col";
           6 "Beam";
           7 "Col";
           8 "Beam";
           9 "Col"];
       
%% Loads       
beams_LL=[1 -30; % Uniformly distributed loads over the beams
          2 -30;
          3 -20;
          4 -20];

elemcols=[];
elembeams=[];
beams=0;
cols=0;
for j=1:nbars
    if type_elem(j,2)=="Beam"
        beams=beams+1;
        elembeams=[elembeams,j];
    elseif type_elem(j,2)=="Col"
        cols=cols+1;
        elemcols=[elemcols,j];
    end
end

% Uniformly distributed loads considering self weight of the elements
qbarxy=zeros(nbars,2);
for i=1:beams
    qbarxy(elembeams(i),2)=1.1*(beams_LL(i,2));
end
   
% Plastic moments of each element's ends
Mp=[7680000 7680000;
    6490000 6490000;
    8363000 8976940;
    5490000 5490000;
    8680000 8680000;
    6490000 6490000;
    9363000 9976940;
    5490000 5490000;
    9680000 9680000]; %Kg-cm

% Initial lateral forces to increment for the Pushover analysis
LatForces=[1500];
            
% Degrees of freedom at which each lateral force is applied (one for
% each force)
dofForces=[4];
Hfloors=[300];

%% Hysteresis curves

%[WDIf,KUDIf,DIf]=hysterCurveEB2DFrames(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
%supports,bc,LatForces,Hfloors,dofForces,0.01,5,1)

%DIf=hysterCurveClough2DFrames(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
%supports,bc,LatForces,Hfloors,dofForces,0.01,8,1)

DIf=hysterCurveTakeda2DFrames(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
supports,bc,LatForces,Hfloors,dofForces,0.01,8,1)