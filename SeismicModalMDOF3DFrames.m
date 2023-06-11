function [fmaxDOF,Mgl,Kgl,T,La,Egv]=SeismicModalMDOF3DFrames...
    (coordxyz,A,unitWeightEl,qbarxyz,eobars,Edof,bc,E,G,J,Iy,Iz,ni,nf,...
    acelxy,g,modal)
% SYNTAX : 
% [fmaxDOF,Mgl,Kgl,T,La,Egv]=SeismicModalMDOF3DFrames...
%  (coordxyz,A,unitWeightEl,qbarxyz,Edof,bc,E,G,J,Iy,Iz,ni,nf,...
%  acelxy,g,modal)
%---------------------------------------------------------------------
%    PURPOSE
%     To compute the global stiffness matrix of a plane frame as well as
%     the global mass matrix to then call the function:
%     "ModalsMDOF2DFrames2" to determine all the structure's modals.
% 
%    INPUT:  coordxyz:          Node coordinates of the structure [x,y,z]
%
%            A:                 Cross-sectional elements' area
%
%            E:                 Modulus of Elasticity of the frame's 
%                               elements
%
%            Iy,Iz:             Cross-sectional inertia with respect to the
%                               local y and z axis of the frame's elements
%
%            bc:                Boundary condition array
%
%            acelxy:            Components of the ground acceleration [x,y]
%                               Only one direction is allowed to be
%                               analysed for each analysis. Therefore, when
%                               analysing with respect to the x global axis
%                               set acelxy = [sax 0], and for the y global
%                               axis acelxy = [0 say],
%
%            modal:             Mode of vibration of interest:
%                               [mode-1,mode-2,...] -> The equivalent 
%                                    inertial forces are computed with the
%                                    contribution of the given modes of 
%                                    vibration
%                               
%                               mode-i -> The equivalent inertial forces are
%                                       computed with the mode of
%                                       vibration inserted
%
%            g:                 gravity acceleration
%
%            unitWeightEl:      unit weight material of each element:
%                               Size: nbars x 1
%
%            qbarxyz:           uniformly distributed loads. 
%                               Size: nbars x 3.
%                               The first column corresponds to the
%                               distributed loads in the local X' direction
%                               the second column to the the loads
%                               distributed in the local Y' direction and
%                               the third column to those in the local Z'
%                               direction
%
%    OUTPUT: La :               Modal of vibration for each DOF. 
%                               Size: Nmodals x 1
%
%            Egv:               DOF's eigenvalues: NDOF x Nmodals
%
%            T :                Structure's periods for each modal      
%
%            fmaxDOF :          Equivalent DOF's forces for the modal
%                               in question
%
%            Mgl:               Global Mass matrix
%            Kgl:               Global Stiffness matrix
%
%--------------------------------------------------------------------

% LAST MODIFIED: L.Verduzco    2023-06-07
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%--------------------------------------------------------------------

nnodes=length(coordxyz(:,1)); nbars=length(E);

%% Stiffness and Mass matrices
Kgl=zeros(6*nnodes);
Mgl=zeros(6*nnodes);
for i=1:nbars      
    ex=[coordxyz(ni(i),1) coordxyz(nf(i),1)];
    ey=[coordxyz(ni(i),2) coordxyz(nf(i),2)];
    ez=[coordxyz(ni(i),3) coordxyz(nf(i),3)];
    ep=[E(i) G(i) A(i) Iy(i) Iz(i) J(i)];
    eo=eobars(i,:);
    
    %% Stiffness matrix
    [Kebar]=beam3e(ex,ey,ez,eo,ep);
    
    %% Mass matrix
    PV=unitWeightEl(i,1)-qbarxyz(i,3)/A(i); % unit weight of each element
                                          % The distributed downward loads
                                          % on the BEAMS are considered.
    [Mebar]=FiniteMassBeams3D(ex,ey,ez,eo,A(i),PV,g,acelxy);
    [Kgl]=assem(Edof(i,:),Kgl,Kebar);
    [Mgl]=assem(Edof(i,:),Mgl,Mebar);
end 

%% Modal analysis
[fmaxDOF,T,La,Egv]=ModalsMDOF3DFrames(Mgl,Kgl,bc,sum(acelxy),modal);

