function [data, elementOut, eType,NodeSets,ElementSets]= abaqusInpRead(fileName)
% read abaqus input file, e.g. fileName = 'myExamInpFile.inp'
% coded by Wan Ji, Wuhan University
% update date 2021/07/15
% last version "readinp", refer to
% https://www.mathworks.com/matlabcentral/fileexchange/95343-read-abaqus-input-file-to-get-the-nodes-and-elements
% something new to this version
% (1) fixed the bug which the function can't read c3d20 element
% (2) offer optional output, that is, the struct output
%%
% If the number of output arguments is only 1, then a structure output is
% activated, the struct includes some important information inside the inp
% file, like "Parts", "Nodes", "Elements", "NodeSets", "ElementSets", "Materials",
% etc.
% Example 01
% fileName = 'myExamInpFile.inp';
% data = abaqusInpRead(fileName);
%%
% If the number of output arguments is 3, then the 1st set of node
% coordinates, the 1st set of elements, and element type will be the output
% data.
% Example 02
% fileName = 'myExamInpFile.inp';
% [node, element, elementType] = abaqusInpRead(fileName);
%% Function body
s = fileread(fileName );
s = lower(s);
s = split(s,'*');
% Substructures of output structure
ElementSets = struct(...
    'ElementSetType',[],...
    'Name',[],...
    'Data',[],...
    'Instance',[],...
    'Generate',[],...
    'Internal',[],...
    'Count', 0);
Parts = struct('Tags', [],...
    'Count', 1,...
    'PartID',[]);
Nodes = struct(...
    'Coordinates', [],...
    'ID', [],...
    'NodeTag',[],...
    'AttachedPartTag',[],...
    'AttachedPartID',[],...
    'Count', 0);
Elements = struct(...
    'ElementType', [],...
    'NodeIDList',[],...
    'ID',[],...
    'ElementTag',[],...
    'AttachedNodeTag',[],...
    'AttachedPartTag', [],...
    'AttachedPartID',[],...
    'Count', 0);
Materials = struct(...
    'MaterialType',[],...
    'MaterialTable', [],...
    'Count', 0);
NodeSets = struct(...
    'NodeSetType',[],...
    'Name',[],...
    'Data',[],...
    'Instance',[],...
    'Generate',[],...
    'Internal',[],...
    'Count', 0);
% the output data
data = struct(...
    'Parts',Parts, ...
    'Nodes', Nodes, ...
    'Elements', Elements, ...
    'Materials',Materials, ...
    'NodeSets', NodeSets, ...
    'ElementSets', ElementSets...
    );
partFlag = true;
for i = 1:1:numel(s)
    [myMatrix, FLS]= getMatrixFromString(s{i}, newline, ',');
    FLS(FLS==char(32)|FLS==char(13)|FLS==newline) = [];
    p = split(FLS,',');
    q = p{1};
    if(isempty(q))
        continue;
    end
    switch lower(q)
        case {'part'}
            if(partFlag)
                partFlag = ~partFlag;
            else
                Parts(1).Count = Parts(1).Count + 1;
            end
            c = Parts(1).Count;
            Parts(c,1).PartID = c;
            partName = split(p{2},'=');
            Parts(c,1).Tags = partName{2};
        case {'node'}
            Nodes(1).Count = Nodes(1).Count + 1;
            c = Nodes(1).Count;
            Nodes(c,1).AttachedPartTag =  Parts(Parts(1).Count).Tags;
            Nodes(c,1).AttachedPartID =  Parts(1).Count;
            idcoor = cell2mat(myMatrix);
            Nodes(c,1).ID = idcoor(:,1);
            Nodes(c,1).Coordinates = idcoor(:,2:end);
            Nodes(c,1).NodeTag = Nodes(1).Count;
        case {'element'}
            Elements(1).Count = Elements(1).Count + 1;
            c = Elements(1).Count;
            Elements(c,1).AttachedPartTag =  Parts(Parts(1).Count).Tags;
            Elements(c,1).AttachedPartID =  Parts(1).Count;
            eType = split(p{2},'=');
            eType = eType{2};
            eNodeNumber = eType(4:end);
            eNodeNumber(double(eNodeNumber)>double('9')|double(eNodeNumber)<double('0')) = [];
            eNodeNumber = str2double(eNodeNumber);
            Elements(c,1).ElementType = eType;
            if(numel(myMatrix)>=2 && (numel(myMatrix{1})+numel(myMatrix{2})-1==eNodeNumber))
                myMat = cell(numel(myMatrix)/2,1);
                for j = 1:1:numel(myMat)
                    myMat{j,1} = [myMatrix{2*j-1,1}, myMatrix{2*j,1}];
                end
            elseif(numel(myMatrix)>=3&& ...
                    (numel(myMatrix{1})+numel(myMatrix{2})+numel(myMatrix{3})-1==eNodeNumber))
                myMat = cell(numel(myMatrix)/3,1);
                for j = 1:1:numel(myMat)
                    myMat{j,1} = [myMatrix{3*j-2,1}, myMatrix{3*j-1,1}, myMatrix{3*j,1}];
                end
            else
                myMat = myMatrix;
            end
            Elements(c,1).NodeIDList = cell2mat(myMat);
            Elements(c,1).ID = Elements(c,1).NodeIDList(:,1);
            Elements(c,1).NodeIDList = Elements(c,1).NodeIDList(:,2:end);
            Elements(c,1).ElementTag = Elements(1).Count;
            Elements(c,1).AttachedNodeTag = Nodes(1).Count;
        case {'nset','ngen'}
            NodeSets(1).Count = NodeSets(1).Count + 1;
            c = NodeSets(1).Count;
            NodeSets(c,1).FirstString = p;
            NodeSets(c,1).Generate = false;
            NodeSets(c,1).Internal = false;
            for j = 1:1:numel(p)
                ep = split(p{j},'=');
                switch ep{1}
                    case {'instance'}
                        NodeSets(c,1).Instance = ep{2};
                    case {'generate'}
                        NodeSets(c,1).Generate = true;
                    case {'internal'}
                        NodeSets(c,1).Internal = true;
                    case {'nset','ngen'}
                        if(numel(ep)>=2)
                            NodeSets(c,1).Name = ep{2};
                        end
                    otherwise
                end
            end
            NodeSets(c,1).Data = cell2mat(myMatrix');
            NodeSets(c,1).SetData = myMatrix;
            NodeSets(c,1).NodeSetType = lower(q);
        case {'elset'}
            ElementSets(1).Count = ElementSets(1).Count + 1;
            c = ElementSets(1).Count;
            ElementSets(c,1).FirstString = p;
            ElementSets(c,1).Generate = false;
            ElementSets(c,1).Internal = false;
            for j = 1:1:numel(p)
                ep = split(p{j},'=');
                switch ep{1}
                    case {'instance'}
                        ElementSets(c,1).Instance = ep{2};
                    case {'generate'}
                        ElementSets(c,1).Generate = true;
                    case {'internal'}
                        ElementSets(c,1).Internal = true;
                    case {'elset'}
                        if(numel(ep)>=2)
                            ElementSets(c,1).Name = ep{2};
                        end
                    otherwise
                end
            end
            ElementSets(c,1).SetData = myMatrix;
            ElementSets(c,1).Data = cell2mat(myMatrix');
            ElementSets(c,1).ElementSetType = lower(q);
        case {'material'}
            Materials(1).Count = Materials(1).Count + 1;
            c = Materials(1).Count;
            Materials(c,1).FirstString = p;
            Materials(c,1).SetData = myMatrix;
            for j = 1:1:numel(p)
                ep = split(p{j},'=');
                switch ep{1}
                    case {'name'}
                        Materials(c,1).Name = ep{2};
                    otherwise
                end
            end
        otherwise
            continue
    end
end
data.Parts = Parts;
data.Nodes = Nodes;
data.Elements = Elements;
data.NodeSets = NodeSets;
data.ElementSets = ElementSets;
data.Materials = Materials;
if(nargout>=2)
    elementOut = data.Elements(1).NodeIDList;
    data = data.Nodes(1).Coordinates;
    eType = Elements(1).ElementType;
    
end
end


function [myMatrix, firstLineString] = getMatrixFromString(s, sepRow, sepColumn)
% get the matrix from a string
ssep = split(s, sepRow);
myMatrix = cell(numel(ssep),1);
myFlag = true(numel(ssep),1);
columnNumber = zeros(numel(ssep),1);
firstLineString = ssep{1};
for i = 1:1:numel(ssep)
    es = ssep{i};
    es(es==char(32)|es==char(13)|es==newline) = [];
    p = split(es, sepColumn);
    myMatrix{i} = zeros(1, numel(p));
    columnNumber(i) = numel(p);
    for j = 1:1:numel(p)
        myMatIJ = str2double(p{j});
        myMatrix{i}(1,j) = myMatIJ;
    end
    myMatrix{i}( isnan(myMatrix{i}))=[];
    if(isempty(myMatrix{i}))
        myFlag(i) = false;
    end
end
myMatrix = myMatrix(myFlag,1);
end