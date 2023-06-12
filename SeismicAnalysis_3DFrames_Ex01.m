% SeismicAnalysis_3DFrames_Ex01
%----------------------------------------------------------------
% PURPOSE 
%    To compute the Static Modal analysis for a 3D Reinforced Concrete
%    Frame.
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-06-10
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------
clc
clear all

nnodes=8;
nbars=8;

%% Materials
% f'c of each element
fpc=[300;
     300;
     300;
     300;
     300;
     300;
     300;
     300];

% Elasticity modulus of each element in function of f'c
E=14000.*sqrt(fpc);

v=[0.2; % Poisson modulus
    0.2;
    0.2;
    0.2;
    0.2;
    0.2;
    0.2;
    0.2];

G=E./(2.*(1+v)); % Shear modulus

%% Geometry/Topology
% cross-section dimensions of each element (rectangular geometry)
dimensions=[30 30;
            25 40;
            30 30;
            25 40;
            25 40;
            30 30;
            25 40;
            30 30];
        
% cross-section area of each element
A=dimensions(:,1).*dimensions(:,2);
        
% Cross-section inertia
Iy=1/12.*dimensions(:,1).*dimensions(:,2).^3;
Iz=1/12.*dimensions(:,2).*dimensions(:,1).^3;

adim=dimensions(:,2).*0.5;
bdim=dimensions(:,1).*0.5;

% Saint Venant constant (polar inertia - Torsion)
J=adim.*bdim.^3.*(16/3-3.36.*bdim./adim.*(1-bdim.^4./(12.*adim.^4)));
      
%% Topology and node coordinates

% OPTION 1: Manually given
% Coordinates of each node
coordxyz=[0 0 0;
          0 0 300;
          0 500 300;
          0 500 0;
          400 0 0;
          400 0 300;
          400 500 300;
          400 500 0];
      
% Connectivity
NiNf=[1 2;
    2 3;
    3 4;
    3 7;
    2 6;
    5 6;
    6 7;
    7 8];

ni=NiNf(:,1);
nf=NiNf(:,2);

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

%% Prescribed boudnary conditions [dof, displacement]
bc=[1 0;
    2 0;
    3 0;
    4 0;
    5 0;
    6 0;
    19 0;
    20 0;
    21 0;
    22 0;
    23 0;
    24 0;
    25 0;
    26 0;
    27 0;
    28 0;
    29 0;
    30 0;
    43 0;
    44 0;
    45 0;
    46 0;
    47 0;
    48 0];

%% Additional data (optional)
type_elem=[1 "Col";
           2 "Beam";
           3 "Col";
           4 "Beam";
           5 "Beam";
           6 "Col";
           7 "Beam";
           8 "Col"];
       
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
       
%% Local z axis of each element
eobars=[0 1 0;
        0 0 1;
        0 1 0;
        0 0 1;
        0 0 1;
        0 1 0;
        0 0 1;
        0 1 0];
    
%% Loads       
beams_LL=[1 -50; % Uniformly distributed loads over the beams
          2 -50;
          3 -50;
          4 -50];

% Assignation of distributed loads on beams (local axis [x',y',z'])

qbarxyz=zeros(nbars,4);
qbarxyz(elembeams',3)=beams_LL(:,2);

%% Mode of vibration of interest
modal=2; % 2 -> acceleration in the x direction
         % 1 -> acceleration in the y direction

%% Seismic response spectrum from the CFE-15

g=981; % gravity acceleration
Fsit=2.4; FRes=3.8; % Factores de sitio y de respuesta
a0_tau=100; % cm/seg^2

ro=0.8; % Redundance factor
alf=0.9; % Irregularity factor
Q=4; % Seismic behaviour factor

Ta=0.1;
Tb=0.6;
Te=0.5; % Structure's period
k=1.5; % Design spectrum slope
Qp=1+(Q-1)*sqrt(Te/(k*Tb)); % Ductility factor

Ro=2.5; % Over-resistance index
R=Ro+1-sqrt(Te/Ta); % Over-resistance factor

sa=-a0_tau*Fsit*FRes/(R*Qp*alf*ro); % Reduced pseudo-acceleration (cm/seg^2)

%% Modal analysis
pvconc=0.0024; % unit weight of concrete
unitWeightElm=zeros(nbars,1)+pvconc;

% Consistent mass method
[fmaxDOF,Mgl,Kgl,T,La,Egv]=SeismicModalMDOF3DFrames...
(coordxyz,A,unitWeightElm,qbarxyz,eobars,Edof,bc,E,G,J,Iy,Iz,NiNf(:,1),...
NiNf(:,2),sa,g,modal);

%% Static structural analysis with seismic forces
np=7; % number of analysis points for the mechanical elements

[edi,eci,displacements,reactions,Ex,Ey,Ez,esbarsnormal,esbarssheary,...
esbarsshearz,esbarstorsion,esbarsmomenty,esbarsmomentz]=StaticLinear3DFrames...
(E,A,Iz,Iy,G,J,bc,fmaxDOF,[1:6*nnodes],NiNf(:,1),NiNf(:,2),...
qbarxyz,5,coordxyz,eobars,1,[],1);

%% Plot of the modal in question and its frequency
Freq=1./T;
zc=0.5*max(coordxyz(:,3));
if length(modal)==1 % If only one modal was entered
    figure(6)
    % Undeformed structure
    grid on
    NoteMode=num2str(modal);
    title(strcat('Eigenmode ','- ',NoteMode))
    elnum=Edof(:,1);
    plotpar=[1,2,1];
    eldraw3(Ex,Ey,Ez,plotpar,elnum)
    
    % Deformed structure
    magnfac=100;
    Edb=extract(Edof,Egv(:,modal));
    plotpar=[1,3,1];
    [magnfac]=eldisp3(Ex,Ey,Ez,Edb,plotpar,magnfac);
    FreqText=num2str(Freq(modal));
    NotaFreq=strcat('Freq(Hz)= ',FreqText);
    text(50,zc,NotaFreq);
end

% ----------------------------- End ----------------------------------