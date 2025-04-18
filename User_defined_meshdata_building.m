%% put your ABAQUS inp model under \meshdata\
%%% Boundary condition can be be pre-defined in inp for processing.
%%%  Only quadrilateral element type is allowed for robust FE computation in optimization process.


clc; clear
cd(fileparts(mfilename('fullpath')));
figure(1)
clf
fileName = [pwd,'\meshdata\Demo-2inpfile.inp']; 
[node, element, elementType,nodesets,elementsets] = abaqusInpRead_revise(fileName);
patch('vertices',node,'faces',element,...
    'facevertexCdata',node(:,1),'facecolor','interp','edgecolor','k');
hold on
colormap jet
axis off
axis equal
%==========BC =================================
% fixednid = [1,2];       
% fixndof =  [1,2,3,4];
% loadnid = [60];
% loadndofs= [120];
% Fvalue = [-1];

%=ABAQUS INP SET id:    1,30,1 --> 1:1:30

%%Demo1
% fixednid =unique([ nodesets(1).Data, nodesets(2).Data ]);
% fixeddofs = sort( [nodesets(1).Data*2-1, nodesets(2).Data*2]) ;
% loadnid = nodesets(3).Data;
% loadndofs = loadnid*2;
% Fvalue=[-2,-1];
% save([pwd,'\meshdata\FEMmeshdata_1.mat'])

%===================================================
%%Demo2
fixednid =unique([ nodesets(1).Data, nodesets(2).Data ]);
fixeddofs = sort( [nodesets(1).Data*2-1, nodesets(2).Data*2]) ;
loadnid = nodesets(3).Data;
loadndofs = loadnid*2-1;
Fvalue= ones( length(loadndofs),1);
save([pwd,'\meshdata\FEMmeshdata_2.mat'])
