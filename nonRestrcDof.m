function [ndof,edof]=nonRestrcDof(nnodes,bc)
%------------------------------------------------------------------------
% Syntax:
% [ndof,edof]=nonRestrcDof(nnodes,bc)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute the DOF vector containing the DOF that are not restricted
%  based on the pre-stablished restricted DOF vector (bc). Through this
%  vector the reduced stiffness matrix could be computed in order to
%  compare the determinants of each structural system at every stage of
%  analysis.
%  
% 
% INPUT:  nnodes:                total number of nodes of the structure
%
%         bc:                    boundary condition vector, containing the
%                                DOF that are restricted in motion or
%                                displacement. Size: [restric-DOF,2] in
%                                format [DOF,displacement]
%
% OUTPUT: edof:                  is vector containing the non-restricted
%                                DOF. Size: [Non-restric-DOF,1]
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-02-18
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

edof_bc=zeros(nnodes*3,1);
edof=[];
for i=1:length(bc(:,1))
    
    edof_bc(bc(i,1))=bc(i,1);
end

for i=1:nnodes*3
    if edof_bc(i)==0
        edof=[edof;i];
    end
end

ndof=length(edof);