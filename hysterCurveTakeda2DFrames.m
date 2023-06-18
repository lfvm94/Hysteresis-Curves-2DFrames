function hysterCurveTakeda2DFrames(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
    supports,bc,seismicForces,Hfloors,dofSeismicForces,dL,ncycles,hcfloor)
%------------------------------------------------------------------------
% Syntax:
% hysterCurveTakeda2DFrames(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
% supports,bc,seismicForces,Hfloors,dofSeismicForces,dL,ncycles,hcfloor)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute the hysteresis curves according to the Takeda's Model 
%  (Load-Deflection diagram of Trilinear degrading stiffness Hysteresis 
%  Model) for a 2D Frame's floor.
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
dispHistFloor1,forceHistFloor1,ncycles);

ky=forceHistFloor1(:,2)./dispHistFloor1(:,2);
uyr=dispHistFloor1(:,2);

%% Pushover to the left
dispHistFloor2=[zerov];
forceHistFloor2=[zerov];

[lambdaLeft,barPlasNodeLeft,dispHistFloor2,forceHistFloor2,supports2,...
plastbars2,mpbar2]=Pushover2DFrames2(qbarxy,A,Mp,E,I,coordxy,ni,nf,...
supports2,bc,-seismicForces,dofSeismicForces,dL,plastbars2,mpbar2,...
dispHistFloor2,forceHistFloor2,ncycles);

uyl=abs(dispHistFloor2(:,2));

%% Unloading points
dg=[];  fg=[];
for i=1:ncycles
    umaxr=max(dispHistFloor1(:,2:i+2)')';
    umaxl=max(abs(dispHistFloor2(:,2:i+2))')';
    krr=ky.*(uyr./umaxr).^0.1;
    krl=ky.*(uyl./umaxl).^0.1;
    
    ulir=dispHistFloor1(:,i+2)-forceHistFloor1(:,i+2)./krr;
    ulil=dispHistFloor2(:,i+2)-forceHistFloor2(:,i+2)./krl;
    if i==1
        % unload and reload points
        dir=[ulir,dispHistFloor2(:,i+1)];
        fir=[zerov,forceHistFloor2(:,i+1)];
        
        % Joining load and unload points
        dg=[dg,dispHistFloor1(:,1:3),dir,dispHistFloor2(:,i+2),ulil];
        fg=[fg,forceHistFloor1(:,1:3),fir,forceHistFloor2(:,i+2),zerov];
    else
        % reload and unload 
        dil=[dispHistFloor2(:,i+1:i+2),ulil];
        fil=[forceHistFloor2(:,i+1:i+2),zerov];
        
        % Joining load and unload points
        dg=[dg,dispHistFloor1(:,i+1:i+2),ulir,dil];
        fg=[fg,forceHistFloor1(:,i+1:i+2),zerov,fil];
    end
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
        ') - Takeda Model'))
    hold on
    grid on
else
    disp('Error. That floor does not exist')
end