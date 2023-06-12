function [Me]=FiniteMassBeams3D(ex,ey,ez,eo,A,Iy,Iz,PV,g)
% SYNTAX : [Me]=FiniteMassBeams3D(ex,ey,ez,eo,A,Iy,Iz,PV,g)
%---------------------------------------------------------------------
%    PURPOSE
%     Compute the mass matrix for a three dimensional beam element. 
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]       element node coordinates
%            ez = [z1 z2]
%
%            eo:                local axis of the element in global
%                               coordinates: 
%                               [component-x, component-y, component-z]
% 
%            A                  Transversal area of element
%            PV                 Volumetric weigth of material element
%            g                  Gravity acceleration

%    OUTPUT: Me : element mass matrix (12 x 12)
%
%--------------------------------------------------------------------
%
% LAST MODIFIED: L.Verduzco    2023-06-10
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%--------------------------------------------------------------------

b=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
L=sqrt(b'*b);  n1=b/L;

lc=sqrt(eo*eo'); n3=eo/lc;

rx=(Iy+Iz)/A;
R=PV*A*L/(420*g);
Mle=R*[140  0   0    0     0   0   70   0   0    0     0      0 ;
        0 156   0    0     0 22*L   0  54   0    0     0  -13*L ;
        0   0  156   0  -22*L  0    0   0   54   0  13*L      0 ;
        0   0   0 140*rx^2 0   0    0   0   0 70*rx^2  0      0 ;
        0   0 -22*L  0  4*L^2  0    0   0 -13*L  0 -3*L^2     0 ;
        0 22*L  0    0     0 4*L^2  0 13*L  0    0     0 -3*L^2 ;
       70   0   0    0     0   0   140  0   0    0     0      0 ;
        0  54   0    0     0 13*L   0 156   0    0     0  -22*L ;
        0   0   54   0  -13*L  0    0   0  156   0  22*L      0 ;
        0   0   0 70*rx^2  0   0    0   0   0 140*rx^2 0      0 ;
        0   0  13*L  0  -3*L^2 0    0   0 22*L   0  4*L^2     0 ;
        0 -13*L 0    0     0 -3*L^2 0 -22*L 0    0     0  4*L^2];
    
n2(1)=n3(2)*n1(3)-n3(3)*n1(2);
n2(2)=-n1(3)*n3(1)+n1(1)*n3(3);
n2(3)=n3(1)*n1(2)-n1(1)*n3(2);

An=[n1';
    n2;
    n3];

G=[  An     zeros(3) zeros(3) zeros(3);
   zeros(3)   An     zeros(3) zeros(3);
   zeros(3) zeros(3)   An     zeros(3);
   zeros(3) zeros(3) zeros(3)   An    ];
  
 Me=G'*Mle*G;