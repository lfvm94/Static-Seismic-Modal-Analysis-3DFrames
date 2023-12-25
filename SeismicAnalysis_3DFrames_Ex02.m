% SeismicAnalysis_3DFrames_Ex02
%----------------------------------------------------------------
% PURPOSE 
%    To compute the Static Modal analysis for a 3D Reinforced Concrete
%    Frame.
%
%----------------------------------------------------------------
%    Notes: There is the option to extract the topology of the 
%           structure from a SAP2000 model by using the SM Toolbox.
%
%           To download the SM Toolbox visit: 
%           https://github.com/RJ-Soft/SM-Toolbox/releases/tag/7.0.2
%----------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-06-01
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------
clc
clear all

%% Topology and node coordinates
%{
% OOPTION 2: Exported from a SAP2000 model (the SM Toolbox is required)
APIDLLPath ='C:\Program Files\Computers and Structures\SAP2000 22\SAP2000v1.dll';
ProgramPath ='C:\Program Files\Computers and Structures\SAP2000 22\SAP2000.exe';

ModelName = 'Frame_Ex01.sdb';
ModelPath = fullfile('C:\Users\luizv\OneDrive\SeismicAnalysis_3DFrames\StaticModalAnalysis\FrameSAP2000_Ex02',ModelName);

[coordxyz,NiNf]=ExtractTopologySAP2000(ProgramPath,APIDLLPath,...
                                            ModelPath);
                                        
coordxyz=coordxyz*2.54; % to change the location coordinates from in to cm
%}
% OPTION 1: Coordinates and connectivity manually given

NiNf=[1,2;2,3;3,4;3,5;5,6;7,8;8,9;9,10;9,11;11,12;8,2;9,3;11,5;2,13;13,14;14,3;14,15;15,5;8,16;16,17;17,9;17,18;18,11;16,13;17,14;18,15;13,19;19,20;20,14;20,21;21,15;16,22;22,23;23,17;23,24;24,18;22,19;23,20;24,21];
coordxyz=[0,0,0; 0,0,118.110236220472; 137.795275590551,0,118.110236220472;137.795275590551,0,0;393.700787401575,0,118.110236220472;393.700787401575,0,0;0,157.480314960630,0;0,157.480314960630,118.110236220472;137.795275590551,157.480314960630,118.110236220472;137.795275590551,157.480314960630,0;393.700787401575,157.480314960630,118.110236220472;393.700787401575,157.480314960630,0;0,0,236.220472440945;137.795275590551,0,236.220472440945;393.700787401575,0,236.220472440945;0,157.480314960630,236.220472440945;137.795275590551,157.480314960630,236.220472440945;393.700787401575,157.480314960630,236.220472440945;0,0,354.330708661417;137.795275590551,0,354.330708661417;393.700787401575,0,354.330708661417;0,157.480314960630,354.330708661417;137.795275590551,157.480314960630,354.330708661417;393.700787401575,157.480314960630,354.330708661417]*2.54;

nnodes=length(coordxyz(:,1));
nbars=39;

% Topology matrix
ni=NiNf(:,1);
nf=NiNf(:,2);
 
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
%% Materials
% f'c of each element
fc1=300;
fpc=zeros(nbars,1)+fc1;

% Elasticity modulus of each element in function of f'c
E=14000.*sqrt(fpc);

v1=0.2;
v=zeros(nbars,1)+v1; % Poisson modulus

G=E./(2.*(1+v)); % Shear modulus

%% Geometry
% cross-section dimensions of each element (rectangular geometry)
dimensions=[30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            25 40;
            25 40;
            25 40;
            30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            25 40;
            25 40;
            25 40;
            30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            30 30;
            25 40;
            30 30;
            25 40;
            30 30;
            25 40;
            25 40;
            25 40];

% cross-section area of each element
A=dimensions(:,1).*dimensions(:,2);

% Cross-section inertia
Iy=1/12.*dimensions(:,1).*dimensions(:,2).^3;
Iz=1/12.*dimensions(:,2).*dimensions(:,1).^3;

adim=dimensions(:,2).*0.5;
bdim=dimensions(:,1).*0.5;

% Saint Venant constant (polar inertia - Torsion)
J=adim.*bdim.^3.*(16/3-3.36.*bdim./adim.*(1-bdim.^4./(12.*adim.^4)));
      
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
    31 0;
    32 0;
    33 0;
    34 0;
    35 0;
    36 0;
    37 0;
    38 0;
    39 0;
    40 0;
    41 0;
    42 0;
    55 0;
    56 0;
    57 0;
    58 0;
    59 0;
    60 0;
    67 0;
    68 0;
    69 0;
    70 0;
    71 0;
    72 0];
       
%% Additional data (optional)
type_elem=[1 "Col";
           2 "Beam";
           3 "Col";
           4 "Beam";
           5 "Col";
           6 "Col";
           7 "Beam";
           8 "Col";
           9 "Beam";
           10 "Col";
           11 "Beam";
           12 "Beam";
           13 "Beam";
           14 "Col";
           15 "Beam";
           16 "Col";
           17 "Beam";
           18 "Col";
           19 "Col";
           20 "Beam";
           21 "Col";
           22 "Beam";
           23 "Col";
           24 "Beam";
           25 "Beam";
           26 "Beam";
           27 "Col";
           28 "Beam";
           29 "Col";
           30 "Beam";
           31 "Col";
           32 "Col";
           33 "Beam";
           34 "Col";
           35 "Beam";
           36 "Col";
           37 "Beam";
           38 "Beam";
           39 "Beam";];
       
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
for i=1:nbars
    if type_elem(i,2)=="Col"
        eobars(i,:)=[0 1 0];
    else
        eobars(i,:)=[0 0 1];
    end
end
    
%% Loads       
beams_LL=-50; % Uniformly distributed loads over the beams

% Assignation of distributed loads on beams
qbarxyz=zeros(nbars,4);
qbarxyz(elembeams',3)=beams_LL;

% Lateral equivalent seismic forces from a modal analysis
nbays=2;

%% Mode of vibration of interest
modal=1; % 2 -> acceleration in the x direction
         % 1 -> acceleration in the y direction

%% Seismic response spectrum from the CFE-15

g=981; % gravity acceleration
DS=1;

%% Modal analysis
pvconc=0.0024; % unit weight of concrete
unitWeightElm=zeros(nbars,1)+pvconc;

% Consistent mass method
[fmaxDOF,Mgl,Kgl,T,La,Egv]=SeismicModalMDOF3DFrames...
(coordxyz,A,unitWeightElm,qbarxyz,eobars,Edof,bc,E,G,J,Iy,Iz,NiNf(:,1),...
NiNf(:,2),DS,g,modal);

%% Static structural analysis with seismic forces
np=7; % number of analysis points for the mechanical elements

[edi,eci,displacements,reactions,Ex,Ey,Ez,esbarsnormal,esbarssheary,...
esbarsshearz,esbarstorsion,esbarsmomenty,esbarsmomentz]=StaticLinear3DFrames...
(E,A,Iz,Iy,G,J,bc,fmaxDOF,[1:6*nnodes],NiNf(:,1),NiNf(:,2),...
qbarxyz,5,coordxyz,eobars,1,[],100);

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