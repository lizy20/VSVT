function Matlab2Abaqus(Nodes,Elements,Elements_Sets,Filename)
fileID = fopen(Filename, 'w');
%Generate Nodes in Input File
fprintf(fileID,'*NODE, NSET=NODE\n');
[NNode, ND]=size(Nodes);
if ND==2  %2D                               
    
for i=1:1:NNode
    
fprintf(fileID,[num2str(i) ', ' num2str(Nodes(i,1)) ', ' num2str(Nodes(i,2)) '\n']); 
end
elseif ND==3  %3D
    
for i=1:1:NNode
    
fprintf(fileID,[num2str(i) ', ' num2str(Nodes(i,1)) ', ' num2str(Nodes(i,2)) ', ' num2str(Nodes(i,3)) '\n']); 
end    
    
end
fprintf(fileID,'\n');
%Generate Elements in Input File
for i=1:1:length(Elements_Sets)
fprintf(fileID,strcat('*ELEMENT, ELSET=',Elements_Sets{i}.Name,', TYPE=',Elements_Sets{i}.Elements_Type,'\n'));
for j=1:1:length(Elements_Sets{i}.Elements) %Loop for the elements in the elements set
IE=Elements_Sets{i}.Elements(j); %Elements indices in elements sets
NNN=[num2str(IE) ', '];
for k=1:1:length(Elements{IE})    
    
NNN=[NNN num2str(Elements{IE}(k)) ', '];
end
NNN=NNN(1:end-2);
fprintf(fileID,[NNN '\n']);
end
fprintf(fileID,'\n');
end
fclose(fileID);
end