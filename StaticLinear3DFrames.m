function [edibars,ecibars,displacements,reactions,Ex,Ey,Ez,esbarsNormal,...
    esbarsShearY,esbarsShearZ,esbarsTorsion,esbarsMomentY,esbarsMomentZ]=...
    StaticLinear3DFrames(E,A,Iz,Iy,G,J,bc,extForce,dofForce,ni,nf,qbarxyzw,...
    np,coordxyz,eobars,plDef,barMech)
% SYNTAX : 
% [edibars,ecibars,displacements,reactions,Ex,Ey,Ez,esbarsNormal,...
%  esbarsShearY,esbarsShearZ,esbarsTorsion,esbarsMomenty,esbarsMomentZ]=...
%  StaticLinear3DFrames(E,A,Iz,Iy,G,J,bc,extForce,dofForce,ni,nf,qbarz,...
%  np,coordxyz,eobars,plDef,barMech)
%---------------------------------------------------------------------
%    PURPOSE
%     To perform a static linear analysis for a 3D Frame.
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
%            J = [polar-inertia;   polar momentum of inertia for all
%                       ...]       elements' cross-sections
% 
%            G = [g_bar;        Shear modulus of elasticity of each element
%               ...]  
%
%            bc:                Boundary condition array
%
%            eobars             local Z axis of each element in the global
%                               system of reference
%
%            qbarxyzw:          uniformly distributed loads. 
%                               Size: nbars x 4.
%                               The first column corresponds to the
%                               distributed loads in the local X' direction
%                               the second column to the the loads
%                               distributed in the local Y' direction,
%                               the third column to those in the local Z'
%                               direction and the fourth corresponds to
%                               torsion distributed loads
%
%            plDef:             parameter that indicates if the plot of the
%                               deformed structure is required or not:
%                               1 -> Plot
%                               else -> Do not plot
%
%            barMech:           list of bars for which it is required to
%                               visualize its mechanical elements
%
%    OUTPUT: 
%
%--------------------------------------------------------------------

% LAST MODIFIED: L.Verduzco    2023-06-10
% Copyright (c)  Faculty of Engineering
%                Autonomous University of Queretaro
%--------------------------------------------------------------------

nbars=length(E);
nnodes=length(coordxyz(:,1));

% External forces
fglobal=zeros(6*nnodes,1);
fglobal(dofForce)=extForce;

% Topology matrix 
Edof=zeros(nbars,13);
for i=1:nbars
    Edof(i,1)=i;
    
    Edof(i,2)=ni(i)*6-5;
    Edof(i,3)=ni(i)*6-4;
    Edof(i,4)=ni(i)*6-3;
    Edof(i,5)=ni(i)*6-2;
    Edof(i,6)=ni(i)*6-1;
    Edof(i,7)=ni(i)*6;
    
    Edof(i,8)=nf(i)*6-5;
    Edof(i,9)=nf(i)*6-4;
    Edof(i,10)=nf(i)*6-3;
    Edof(i,11)=nf(i)*6-2;
    Edof(i,12)=nf(i)*6-1;
    Edof(i,13)=nf(i)*6;
end

Kglobal=zeros(6*nnodes);
collecKebars=zeros(12*nbars,12);
 for i=1:nbars      
     
     ex=[coordxyz(ni(i),1) coordxyz(nf(i),1)];
     ey=[coordxyz(ni(i),2) coordxyz(nf(i),2)];
     ez=[coordxyz(ni(i),3) coordxyz(nf(i),3)];
     
     eo=eobars(i,:);
     ep=[E(i) G(i) A(i) Iy(i) Iz(i) J(i)];
     eq=qbarxyzw(i,:);
     
     [Kebar,febar]=beam3e(ex,ey,ez,eo,ep,eq);
     collecKebars((i-1)*12+1:i*12,:)=Kebar;
     
     [Kglobal,fglobal]=assem(Edof(i,:),Kglobal,Kebar,fglobal,febar);
     
 end     
[displacements,reactions]=solveq(Kglobal,fglobal,bc);
Ed=extract(Edof,displacements);

ex=coordxyz(:,1);
ey=coordxyz(:,2);
ez=coordxyz(:,3);

Ex=zeros(nbars,2);
Ey=zeros(nbars,2);
Ez=zeros(nbars,2);
for j=1:nbars
    Ex(j,1)=ex(Edof(j,7)/6);
    Ex(j,2)=ex(Edof(j,13)/6);

    Ey(j,1)=ey(Edof(j,7)/6);
    Ey(j,2)=ey(Edof(j,13)/6);
    
    Ez(j,1)=ez(Edof(j,7)/6);
    Ez(j,2)=ez(Edof(j,13)/6);
end
 
edibars=zeros(np*nbars,4);
 for i=1:nbars
     ep=[E(i) G(i) A(i) Iy(i) Iz(i) J(i)];
     eq=qbarxyzw(i,:);
     
     [es_bar,edi,eci]=beam3s(Ex(i,:),Ey(i,:),Ez(i,:),eobars(i,:),ep,...
         Ed(i,:),eq,np);
     esbarsNormal(:,i)=es_bar(:,1);
     esbarsShearY(:,i)=es_bar(:,2);
     esbarsShearZ(:,i)=es_bar(:,3);
     esbarsTorsion(:,i)=es_bar(:,4);
     esbarsMomentY(:,i)=es_bar(:,5);
     esbarsMomentZ(:,i)=es_bar(:,6);
     
     ecibars(:,i)=eci;
     edibars((i-1)*np+1:(i*np),:)=edi;
 end
 

%% Mechanical element diagrams

if isempty(barMech)==0
    %% TORSION - MOM Y - MOM Z
    figure(1)
    plot(ecibars(:,barMech),esbarsTorsion(:,barMech),'b -');
    hold on
    ylabel('Magnitude')
    title('Mechanical elements')
    xlabel('X')
    legend('Torsion')
    
    figure(1)
    plot(ecibars(:,barMech),esbarsMomentY(:,barMech),'k -',...
        'linewidth',0.1,'MarkerFaceColor','black','DisplayName','My')
    hold on

    figure(1)
    plot(ecibars(:,barMech),esbarsMomentZ(:,barMech),'g -',...
        'linewidth',0.1,'MarkerFaceColor','green','DisplayName','Mz')
    hold on

    %% NORMAL - CORTANTE Y - CORTANTE Z
    plot(ecibars(:,barMech),esbarsShearZ(:,barMech),'b -');
    hold on
    ylabel('Magnitude')
    title('Mechanical elements')
    xlabel('X')
    legend('Vz')

    figure(2)
    plot(ecibars(:,barMech),esbarsShearY(:,barMech),'r -',...
        'linewidth',0.1,'MarkerFaceColor','red','DisplayName','Vy')
    hold on

    figure(2)
    plot(ecibars(:,barMech),esbarsNormal(:,barMech),'g -',...
        'linewidth',0.1,'MarkerFaceColor','green','DisplayName','Normal')
    hold on
end

%% Deformed-Undeformed structure
if plDef==1
    figure(3)
    % Undeformed structure
    plotpar=[1,2,1];
    elnum=Edof(:,1);
    eldraw3(Ex,Ey,Ez,plotpar,elnum)
    
    % Deformed structure
    magnfac=100;
    Ed=extract(Edof,displacements);
    plotpar=[1,3,1];
    [magnfac]=eldisp3(Ex,Ey,Ez,Ed,plotpar,magnfac);
    
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(strcat('Deformed-undeformed structure - Positive Forces; ',...)
    ' Scale: 1- ', num2str(magnfac)));
end