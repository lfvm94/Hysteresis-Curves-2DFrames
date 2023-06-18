function [historyIncLoad,barPlasNode,dispHistFloor,forceHistFloor,support,...
plastbars,mpbar]=Pushover2DFrames2(qbarxy,A,Mp,E,I,coordxy,ni,nf,support,...
bc,seismicforces,dofForces,dload,plastbars,mpbar,dispHistFloor,...
forceHistFloor,ncycles)

%------------------------------------------------------------------------
% [historyIncLoad,pdriftDI,driftDI,defBasedDI,maxDisplacement,...
%    barPlasNode,dispHistFloor]=Pushover2DFrames2(qbary,A,Mp,E,I,...
%    coordxy,ni,nf,support,bc,seismicforces,dofForces,dload,ncycles)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute a static non-linear pushover analysis of a plane frame
%  
% 
% INPUT:  A = [area_bar;
%               ...]                 area of all elements
%
%         Mp = [Mpi Mpj;             Plastic Moment for each member 
%               ... ]                (i) initial node, (j) final node
%
%         E = [e_bar;                Elasticity modulus of each element
%               ...]                    
%
%         I = [inertia_bar;         in-plane inertia for all elements'
%                       ...]         cross-section
%                                    
%         coordxy = [coordx coordy;  node coordinates for all nodes
%                       ...];
%
%         ni                         list of initial nodes of all bars,
%         nf                         list of final nodes of all bars:
%                                         size = [nbars,1]
% 
%         qbarxy = [bar, load;
%                   ..    .. ]       uniform distributed load acting 
%                                    downwards (only for beams)
%
%         support = [i, j]           support at each bar's end
%                                    options: "Art" or "Fixed"
%                                    (i) initial node, (j) final node
%
%         bc                         restricted dof
%
%         seismicForces = [f(1);]    lateral forces per floor:
%                          f(n);]    size = [nfloors,1]
%
%         Hfloor = [h(1);            Height of each floor from bottom
%                    h(n)]           to top: size = [nfloors,1]
%
%         dofForces = [dof-f(1),     dof at which the lateral forces are
%                       dof-f(n)]    applied - global
%
% OUTPUT: historyIncLoad             history of incremental load factors at
%                                    at which plastic moments are reached
%
%         pdriftDI                   Plastic inter-story drift Damage 
%                                    Index per floor: size = [nfloors,1]
%
%         driftDI                    Inter-story drift Damage 
%                                    Index per floor: size = [nfloors,1]
%
%         defBasedDI                 Deformation based Damage Index
%                                    per floor size = [nfloors,1]
%
%         maxDisplacement            Max absoloute lateral displacement
%                                    for each floor: size = [nfloors,1]
% 
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-01-18
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

nbars=length(E);
nnodes=length(coordxy(:,1));

% Topology matrix 
Edof=zeros(nbars,7);
for i=1:nbars
    Edof(i,1)=i;
    Edof(i,2)=ni(i)*3-2;
    Edof(i,3)=ni(i)*3-1;
    Edof(i,4)=ni(i)*3;
    
    Edof(i,5)=nf(i)*3-2;
    Edof(i,6)=nf(i)*3-1;
    Edof(i,7)=nf(i)*3;
end
    
[ndof,edof]=nonRestrcDof(nnodes,bc);

barPlasNode=[];
iter_collection=[];
historyIncLoad=[];
looping=0;
incLoad=1.0;
iteration=0;
nplas=0;
while looping==0
    Kglobal=zeros(3*nnodes);
    
    fglobal=zeros(3*nnodes,1);
    fglobal(dofForces)=seismicforces*incLoad;
    
    elmMat=zeros(6*nbars,6);

    Ex=zeros(nbars,2);
    Ey=zeros(nbars,2);
    
    for i=1:nbars 
        ex=[coordxy(ni(i),1) coordxy(nf(i),1)];
        ey=[coordxy(ni(i),2) coordxy(nf(i),2)];
        
        Ex(i,:)=ex;
        Ey(i,:)=ey;
        
        ep=[E(i) A(i) I(i)];
         
        eq=[0 qbarxy(i,2)];

        [Kebar,febar]=beam2e(ex,ey,ep,eq); % This is a CALFEM function
                                           % Download at: 
                                           % https://www.byggmek.lth.se/english/calfem/
        
        if support(i,2)=="Fixed" && support(i,3)=="Art"
            Mpl=mpbar(i,2);
            eq=[0 qbarxy(i,2) Mpl];
            [Kebar,febar]=beamArt2e(ex,ey,ep,eq,1);
             
         elseif support(i,2)=="Art" && support(i,3)=="Fixed"

             Mpl=mpbar(i,1);
             eq=[0 qbarxy(i,2) Mpl];
             [Kebar,febar]=beamArt2e(ex,ey,ep,eq,2);
             
        elseif support(i,2)=="Art" && support(i,3)=="Art"

             Mp1=mpbar(i,1);
             Mp2=mpbar(i,2);
             eq=[0 qbarxy(i,2) [Mp1,Mp2]];
             [Kebar,febar]=beamArt2e(ex,ey,ep,eq,3);
        end
        fbe(:,i)=febar; % storing elemental forces for further use
        
        elmMat((i-1)*6+1:6*i,:)=Kebar; % storing Ke of bars for further use
        
        % Assembling global stiffness matrix
        [Kglobal,fglobal]=assem(Edof(i,:),Kglobal,Kebar,fglobal,febar);
     end     

    globalKreduced=Kglobal(edof,edof);
    if iteration==0
        det_Kred=det(globalKreduced);
    else
        det_Kred_post=det(globalKreduced);
    end
    % Solving the system of equations
    [Uglobal,Reactions]=solveq(Kglobal,fglobal,bc);
   
    % --- computation of mechanic elements at the ends of bars --- %
    Ed=extract(Edof,Uglobal);
    for i=1:nbars
        
        es_bar=beam2s(Ex(i,:),Ey(i,:),[E(i) A(i) I(i)],Ed(i,:),...
            [0 qbarxy(i,2)],2);
        
        ue=Uglobal(Edof(i,2:7));
        ke=elmMat((i-1)*6+1:6*i,:);

        fe=-(ke*ue-fbe(:,i));
        
        reac_bars(:,i)=fe;
    end

    iteration=iteration+1;
    current_plas=0; % to register if there is a plastification in the 
                    % current run

    plastified_bars=zeros(nbars,1); % to register which bars are plastified
                                    % in the current run (if any)
                                    
    for i=1:nbars
        % Detect if any end of this bar (i) has been plastified
        bar_plas_check=0;
        if plastbars(1,i)~=0 || plastbars(2,i)~=0
            bar_plas_check=bar_plas_check+1;
        end

        if bar_plas_check==1
            % Detect if the other end has been also plastified
            if abs(reac_bars(3,i))>=Mp(i,1) && ...
               abs(reac_bars(6,i))>=Mp(i,2) % if both ends are plastified
           
                if plastbars(1,i)==1 && plastbars(2,i)==0 
                    % The bar is currently Art-Fixed and will be Art-Art
                    current_plas=1;

                    plastbars(2,i)=1;

                    % Change condition Art-Fixed to Art-Art
                    support(i,3)="Art";

                    % Equivalent plastic moments
                    mplas=reac_bars(6,i);
                    mpbar(i,2)=mplas;
    
                    plastified_bars(i,1)=2;
                
                elseif plastbars(1,i)==0 && plastbars(2,i)==1
                    % The bar is currently Fixed-Art and will be Art-Art
                    current_plas=1;
                    plastbars(1,i)=1;
                    
                    % Change condition Fixed-Art to Art-Art
                    support(i,2)="Art";

                    % Equivalent plastic moments
                    mplas=reac_bars(3,i);
                    mpbar(i,1)=mplas;
                    
                    plastified_bars(i,1)=1;
                end
            end
        elseif bar_plas_check==0

            if abs(reac_bars(3,i))>=Mp(i,1) && ...
                    abs(reac_bars(6,i))<Mp(i,2)

                current_plas=1;
                mplas=reac_bars(3,i);
                plastbars(1,i)=1;

                % change condition to Art
                support(i,2)="Art";

                % Equivalent plastic moments
                mpbar(i,1)=mplas;
                
                plastified_bars(i,1)=1;

            elseif abs(reac_bars(6,i))>=Mp(i,2) && ...
                    abs(reac_bars(3,i))<Mp(i,1)
                current_plas=1;
                plastbars(2,i)=1;
                mplas=reac_bars(6,i);

                % change condition to Fixed-Art
                support(i,3)="Art";

                % Equivalent plastic moments
                mpbar(i,2)=mplas;
                
                plastified_bars(i,1)=2;
                
            elseif abs(reac_bars(6,i))>=Mp(i,2) && ...
                    abs(reac_bars(3,i))>=Mp(i,1)
                
                current_plas=1;
                plastbars(2,i)=1; 
                plastbars(1,i)=1;
                mplas1=reac_bars(3,i); % storing plastic moments
                mplas2=reac_bars(6,i); % for the next iteration
                
                % change condition to Art-Art
                support(i,2)="Art";
                support(i,3)="Art";

                % Equivalent plastic moments
                mpbar(i,1)=mplas1;
                mpbar(i,2)=mplas2;
                                
                % To have register that both element's ends were
                % plastified in flexure
                plastified_bars(i,1)=3;

            end
        end
    end
    if current_plas==0
        % Updating loads for the next iteration (in case there is one)
        incLoad=incLoad+dload;
    else
        nplas=nplas+1;
        historyIncLoad=[historyIncLoad,incLoad];
        iter_collection=[iter_collection,iteration];
                
        barPlasNode=[barPlasNode,plastified_bars];
        
        disp_iter=Uglobal(dofForces);
        force_iter=seismicforces*incLoad;
        
        dispHistFloor=[dispHistFloor,disp_iter];                
        forceHistFloor=[forceHistFloor,force_iter];
        
        % Updating loads for the next iteration (in case there is one)
        incLoad=incLoad+dload;
        % incLoad=1;
    end
    if iteration>1 % to verify is a next run is required
        if nplas==ncycles+1
            break;
        else
            continue;
        end
    end

end

%---------------------------------end----------------------------------