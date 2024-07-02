function [WDIf,KUDIf,DIf]=hysterCurveEB2DFrames(qbarxy,A,Mp,E,I,...
    coordxy,ni,nf,supports,bc,seismicForces,Hfloors,dofSeismicForces,dL,...
    ncycles,hcfloor)

%------------------------------------------------------------------------
% Syntax:
% [WDIf,KUDIf,DIf]=hysterCurveEB2DFrames(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
% supports,bc,seismicForces,Hfloors,dofSeismicForces,dL,ncycles,hcfloor)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute the hysteresis curves according to the Energy Based Model 
%  of Sucuoglu and Erberik (Load-Deflection diagram of Energy-Based)
%  Hysteresis Model) for a 2D Frame's floor.
%  
% 
% INPUT:  ni,nf:                 are the vectors containing the initial
%                                and final nodes for each element. Size: 
%                                [nbars,1] for each
%
%         qbarxy:                are the distributed loads on the
%                                elements. Size: nbars x 2. The 1st column 
%                                is for the loads in the local x
%                                direction and 2nd column for the loads in
%                                the local y direction of each element.
%
%         Mp = [Mpi Mpj;         Plastic Moment for each member 
%               ... ]            (i) initial node, (j) final node
%
%         E:                     is the vector containing the elasticity
%                                modulus for each bar
%
%         A:                     is the vector containing the cross-section
%                                area of each element
%
%         I:                     is the vector containing the cross-section
%                                inertia momentum of each element
%
%         bc:                    is the array containing the boundary 
%                                conditions for the respective prescribed 
%                                (or restricted) DOF. Size=[nRestrictedDOF,2]
%                                in format [DOF,prescribed-displacement]
%
%         coordxy:               is the array containing the node coordinates.
%                                Size = [nNodes,2] in format [xi,yi]
%
%         support = [i, j]       support at each bar's end
%                                options: "Art" or "Fixed"
%                                (i) initial node, (j) final node
%
%         seismicForces=[f(1);]  lateral forces per floor:
%                        f(n);]  size = [nfloors,1]
%
%         Hfloor = [h(1);        Height of each floor from bottom
%                   h(n)]        to top: size = [nfloors,1]
%
%         dofForces=[dof-f(1),   dof at which the lateral forces are
%                    dof-f(n)]   applied - global
%
%         dL:                    incremental load step of analysis
%
%         ncycles:               number of cycles

%         hcfloor:               floor of interest for the computation of
%                                the hysteresis curves
%
% OUTPUT: WDIf:                  Wang's DI. It is based on the original 
%                                Park-Ang DI (see Doc)
%
%         KUDIf                  Kumar & Usami's DI. It is also a further
%                                version of the original Park-Ang DI (see
%                                Doc)
%
%         DIf:                   is the Low-Fatigue Damage Index for the
%                                floor in question (hcfloor)
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-06-17
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

nfloors=length(Hfloors);
nbars=length(E);
plastbars=zeros(2,nbars);
mpbar=zeros(nbars,2); % to save the plastic moments at each articulation
                     % of each bar as a plastification occurs

plastbars2=zeros(2,nbars);
mpbar2=zeros(nbars,2); % to save the plastic moments at each articulation
                     % of each bar as a plastification occurs
supports2=supports;
%% Pushover to the right
zerov=zeros(nfloors,1);
dispHistFloor1=[zerov];
forceHistFloor1=[zerov];

[lambdaRight,barPlasNodeRight,dispHistFloor1,forceHistFloor1,supports,...
plastbars,mpbar]=Pushover2DFrames2(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
supports,bc,seismicForces,dofSeismicForces,dL,plastbars,mpbar,...
dispHistFloor1,forceHistFloor1,2*ncycles-1);

kyr=forceHistFloor1(:,2)./dispHistFloor1(:,2);

%% Pushover to the left
dispHistFloor2=[zerov];
forceHistFloor2=[zerov];

[lambdaLeft,barPlasNodeLeft,dispHistFloor2,forceHistFloor2,supports2,...
plastbars2,mpbar2]=Pushover2DFrames2(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
supports2,bc,-seismicForces,dofSeismicForces,dL,plastbars2,mpbar2,...
dispHistFloor2,forceHistFloor2,2*ncycles-1);

kyl=forceHistFloor2(:,2)./dispHistFloor2(:,2);

%% Unloading points
dg=[];  fg=[];
for i=1:ncycles
    ulir=dispHistFloor1(:,i+2)-forceHistFloor1(:,i+2)./kyr;
    ulil=dispHistFloor2(:,i+2)-forceHistFloor2(:,i+2)./kyl;
    if i==1
        dir=[ulir,dispHistFloor2(:,i+1)];
        fir=[zerov,forceHistFloor2(:,i+1)];
        
        dg=[dg,dispHistFloor1(:,1:3),dir,dispHistFloor2(:,i+2),ulil];
        fg=[fg,forceHistFloor1(:,1:3),fir,forceHistFloor2(:,i+2),zerov];
    else
        dil=[dispHistFloor2(:,i+2:i+3),ulil];
        fil=[forceHistFloor2(:,i+2:i+3),zerov];
        
        dg=[dg,dispHistFloor1(:,2*i:2*i+1),ulir,dil];
        fg=[fg,forceHistFloor1(:,2*i:2*i+1),zerov,fil];
    end
end

%% Cumulative Damage Indices
Hf=sum(Hfloors(1:hcfloor)); % accumulated height of the floor in question
beta=0.1; % Factor that considers the cyclic loading effect
du=0.04*Hf; 

% For the right half-cycle:
dy=dispHistFloor1(hcfloor,2); 
Fy=forceHistFloor1(hcfloor,2);
dm=max(dispHistFloor1(hcfloor,:)); 
c=1;
%% Wang Park-Ang DI
sumE=0;
Eiv=[];
Et=2*Fy*dy;
dim=dispHistFloor1(hcfloor,3);
for i=1:ncycles
    if i==1
        Ei=2.5*Fy*(dim-dy);
    else
        ki=forceHistFloor1(hcfloor,2*i)/...
            (dispHistFloor1(hcfloor,2*i)-dg(hcfloor,6*i-5));
        Fn=ki*dim;
        Ei=2.5*Fn*(dim-dy);
    end
    Eiv=[Eiv,Ei];
    sumEj=sum(Eiv);
    
    betai=(Ei/(Et-sumEj))^c;
    sumE=sumE+betai*Ei;
end
WDIf=(1-beta)*(dm-dy)/(du-dy)+beta*sumE/(Fy*(du-dy));

%% Kumar & Usami DI
sumN1=0; sumNh=0;
for i=1:ncycles
    sumN1=sumN1+((dm-dy)/(du-dy))^c;
    if i==1
        Ei=2.5*Fy*(dim-dy);
    else
        ki=forceHistFloor1(hcfloor,2*i)/...
            (dispHistFloor1(hcfloor,2*i)-dg(hcfloor,6*i-5));
        Fn=ki*dim;
        Ei=2.5*Fn*(dim-dy);
    end
    sumNh=sumNh+(Ei/(Fy*(du-dy)))^c;
end
KUDIf=(1-beta)*sumN1+beta*sumE/(Fy*(du-dy));

%% Low-Cycle Fatigue Damage Index
DIf=0;
duy=du-dy;
dmy=dm-dy;
b=1.6;
for i=1:ncycles
    DIf=DIf+(dmy/duy)^b;
end

%% Plots
if hcfloor<=nfloors && hcfloor>0
    % Ductility curves
    floorText(1,:)=strcat('Floor ',num2str(hcfloor));
    figure(3)
    plot(dg(hcfloor,:),...
         fg(hcfloor,:),'k -','LineWidth',1.8)
    legend(floorText(1,:))
    
    % Labels
    xlabel('Lateral relative displacement')
    ylabel('Load')
    title(strcat('Hysteresis Curves for Floor (',num2str(hcfloor),...
        ') - Energy Based Model'))
    hold on
    grid on
else
    disp('Error. That floor does not exist')
end