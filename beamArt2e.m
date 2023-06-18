function [Ke,fe]=beamArt2e(ex,ey,ep,eq,art);
% Ke=beamArt2e(ex,ey,ep)
% [Ke,fe]=beamArt2e(ex,ey,ep,eq)
%---------------------------------------------------------------------
%    PURPOSE
%     Compute the stiffness matrix for a two dimensional beam element. 
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]       element node coordinates
%
%            ep = [E A I]       element properties
%                                  E: Young's modulus
%                                  A: Cross section area
%                                  I: Moment of inertia
%
%            eq = [qx qy Mp]    distributed loads, local directions
% 
%            art = 1, 2, 3      articulation condition at the ends
%                                  1. Fixed - Articulated
%                                  2. Articulated - Fixed
%                                  3. Articulated - Articulated
%
%    OUTPUT: Ke : element stiffness matrix (6 x 6)
%
%            fe : element load vector (6 x 1)
%--------------------------------------------------------------------

% LAST MODIFIED: L.F.Verduzco    2023-05-31
% Copyright (c)  Faculty of Engineering.
%                Autonomous University of Queretaro
%-------------------------------------------------------------
b=[ ex(2)-ex(1); ey(2)-ey(1) ];
L=sqrt(b'*b);  n=b/L;

E=ep(1);  A=ep(2);  I=ep(3);

qx=0; qy=0; Mp=0; if nargin>3; qx=eq(1); qy=eq(2); Mp=eq(3); end

if art==1
    Kle=[E*A/L   0            0      -E*A/L      0          0 ;
         0   3*E*I/L^3   3*E*I/L^2   0   -3*E*I/L^3         0 ;
         0   3*E*I/L^2    3*E*I/L    0   -3*E*I/L^2         0 ;
       -E*A/L    0            0      E*A/L       0          0 ;
         0   -3*E*I/L^3 -3*E*I/L^2   0   3*E*I/L^3          0 ;
         0       0            0      0           0          0];

    fle=L*[qx/2 5/8*qy qy*L/8 qx/2 3/8*qy 0]'+...
    	Mp*[0 -3/(2*L) 1/2 0 3/(2*L) 1]';
    
elseif art==2

    Kle=[E*A/L   0            0   -E*A/L         0          0 ;
         0   3*E*I/L^3        0      0   -3*E*I/L^3  3*E*I/L^2;
         0       0            0      0   	     0          0 ;
       -E*A/L    0            0    E*A/L         0          0 ;
         0   -3*E*I/L^3 	  0      0   3*E*I/L^3  -3*E*I/L^2;
         0   3*E*I/L^2    	  0      0   -3*E*I/L^2   3*E*I/L];

    fle=L*[qx/2 3/8*qy 0 qx/2 5/8*qy -qy*L/8]'+...
    	Mp*[0 3/(2*L) 1 0 -3/(2*L) 1/2]';
    
elseif art==3
    Mp1=eq(3); Mp2=eq(4);
    
    Kle=[E*A/L   0            0      -E*A/L      0          0 ;
         0       0            0      0           0          0 ;
         0       0            0      0   	     0          0 ;
       -E*A/L    0            0      E*A/L       0          0 ;
         0       0 	          0      0           0          0 ;
         0       0    	      0      0           0          0];
    
    fle=L*[qx/2 qy/2 0 qx/2 qy/2 0]'+...
    	[0 (Mp1+Mp2)/L Mp1 0 -(Mp1+Mp2)/L Mp2]';
end

G=[n(1) n(2)  0    0    0   0;
-n(2) n(1)  0    0    0   0;
  0    0    1    0    0   0;
  0    0    0   n(1) n(2) 0;
  0    0    0  -n(2) n(1) 0;
  0    0    0    0    0   1];

Ke=G'*Kle*G;   fe=G'*fle; 
%--------------------------end--------------------------------