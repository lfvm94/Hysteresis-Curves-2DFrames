% Example_HysterCurves_Frame_02
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

nnodes=15;
nbars=20;

%% Materials
% f'c of each element
fpc=[300;
     300;
     300;
     300;
     280;
     280;
     280;
     280;
     300;
     300;
     300;
     300;
     280;
     280;
     280;
     280;
     300;
     300;
     300;
     300];

% Elasticity modulus of each element in function of f'c
for i=1:nbars
    E(i)=14000*(fpc(i))^0.5;
end

%% Geometry/Topology
dimensions=[35 70;
            65 70;
            55 70;
            55 55;
            55 110;
            55 110;
            35 70;
            35 70;
            65 105;
            65 90;
            65 90;
            65 205;
            35 70;
            35 70;
            35 70;
            35 70;
            65 70;
            55 70;
            45 60;
            35 60];

for i=1:nbars
    A(i)=dimensions(i,1)*dimensions(i,2);
    I(i)=1/12*dimensions(i,1)*dimensions(i,2)^3;
end

% coordinates of each node for each bar
coordxy=[0 -150;
         0 400;
         0 800;
         0 1200;
         0 1600;
         600 1600;
         600 1200;
         600 800;
         600 400;
         600 -150;
         1200 -150;
         1200 400;
         1200 800;
         1200 1200;
         1200 1600]; 
        
%%%---- Initial-final node of each bar -----%%%

ni=[1;2;3;4;5;4;3;2;10;9;8;7;6; 7; 8; 9; 11;12;13;14];
nf=[2;3;4;5;6;7;8;9;9; 8;7;6;15;14;13;12;12;13;14;15];

L=sqrt((coordxy(nf,1)-coordxy(ni,1)).^2+...
      (coordxy(nf,2)-coordxy(ni,2)).^2); % bar-length vector

% prescribed boudnary conditions [DOF, displacement]
bc=[1 0;
    2 0;
    3 0;
    28 0;
    29 0;
    30 0;
    31 0;
    32 0;
    33 0];

%% Loads  
type_elem=[1 "Col";
           2 "Col";
           3 "Col";
           4 "Col";
           5 "Beam";
           6 "Beam";
           7 "Beam";
           8 "Beam";
           9 "Col";
           10 "Col";
           11 "Col";
           12 "Col";
           13 "Beam";
           14 "Beam";
           15 "Beam";
           16 "Beam";
           17 "Col";
           18 "Col";
           19 "Col";
           20 "Col"];
       
beams_LL=[1 -80; % Uniformly distributed loads over the beams
          2 -80;
          3 -80;
          4 -80;
          5 -80;
          6 -80;
          7 -80;
          8 -80];

elem_cols=[];
elem_beams=[];
beams=0;
cols=0;
for j=1:nbars
    if type_elem(j,2)=="Beam"
        beams=beams+1;
        elem_beams=[elem_beams,j];
    elseif type_elem(j,2)=="Col"
        cols=cols+1;
        elem_cols=[elem_cols,j];
    end
end


supports=[1 "Empotrado" "Empotrado";
       2 "Empotrado" "Empotrado";
       3 "Empotrado" "Empotrado";
       4 "Empotrado" "Empotrado";
       5 "Empotrado" "Empotrado";
       6 "Empotrado" "Empotrado";
       7 "Empotrado" "Empotrado";
       8 "Empotrado" "Empotrado";
       9 "Empotrado" "Empotrado";
       10 "Empotrado" "Empotrado";
       11 "Empotrado" "Empotrado";
       12 "Empotrado" "Empotrado";
       13 "Empotrado" "Empotrado";
       14 "Empotrado" "Empotrado";
       15 "Empotrado" "Empotrado";
       16 "Empotrado" "Empotrado";
       17 "Empotrado" "Empotrado";
       18 "Empotrado" "Empotrado";
       19 "Empotrado" "Empotrado";
       20 "Empotrado" "Empotrado"];

% Uniformly distributed loads considering self weight of the elements
qbarxy=zeros(nbars,2);
for i=1:beams
    qbarxy(elem_beams(i),2)=1.1*(beams_LL(i,2));
    
end
    
% Plastic moments of each element's ends
Mp=[9680000 9680000;
    8490000 8490000;
    8363000 8976940;
    7490000 7490000;
    5680000 5680000;
    7363000 7976940;
    8363000 8976940;
    9490000 9490000;
    12680000 12680000;
    11490000 11490000;
    10363000 10976940;
    9490000 9490000;
    5680000 5680000;
    7363000 7976940;
    8363000 8976940;
    9490000 9490000;
    9680000 9680000;
    8490000 8490000;
    8363000 8976940;
    7490000 7490000]; %Kg-cm

% Initial lateral forces to increment for the Pushover analysis
LatForces=[1500;
            2000;
            2500;
            3000]; 

% Degrees of freedom at which each lateral force is applied (one for
% each force)
dofForces=[4;7;10;13];
Hfloors=[400;400;400;400];

%% Hysteresis curves
[WDIf,KUDIf]=hysterCurveEB2DFrames(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
supports,bc,LatForces,Hfloors,dofForces,0.01,3,1)

%hysterCurveClough2DFrames(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
%supports,bc,LatForces,Hfloors,dofForces,0.01,4,4)

%hysterCurveTakeda2DFrames(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
%supports,bc,LatForces,Hfloors,dofForces,0.01,5,4)

