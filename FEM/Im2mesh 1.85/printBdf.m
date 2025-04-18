function printBdf( nodecoor_list, ele_cell, precision_nodecoor )
% printBdf: print the nodes and elements into Inp file 'test.bdf'
%
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, Oct 2020
% Cite As
%   Jiexian Ma (2020). Im2mesh (2D image to triangular meshes) (https://ww
%   w.mathworks.com/matlabcentral/fileexchange/71772-im2mesh-2d-image-to-t
%   riangular-meshes), MATLAB Central File Exchange. Retrieved May 21, 202
%   0.
%   
    
    num_node = size( nodecoor_list, 1 );
    num_phase = length( ele_cell );
    
    % ------------------------------------------------------------------------
    % format of number
    % field width of node numbering
    width_node_num = 1 + floor( log10( num_node ) );         % 18964 -> 5
    if width_node_num > 16
        error('more than 16 digits')
    end
    
    num_digits_of_int_part = 1 + floor( log10( max(nodecoor_list(end,2:3)) ) );
                                                             % 182.9 -> 3
    if num_digits_of_int_part + precision_nodecoor + 1 > 16
        error('more than 16 digits')
    end
    format_node_coor = [ '%.', num2str( precision_nodecoor ), 'f' ];
                                                             % '%.(precision)f'
    
    % ------------------------------------------------------------------------
	fid=fopen('test.bdf','wW');
    % ------------------------------------------------------------------------
    fprintf( fid, 'BEGIN BULK\n');
    
    % print node
    % GRID*,3,,0.5000,2.5000
    fprintf( fid, ...
            [ 'GRID*,%d,,', format_node_coor, ',', format_node_coor, '\n'], ...
            nodecoor_list' );

    % ------------------------------------------------------------------------
    % renumber global element numbering in ele_cell{i}(:,1)
    count = 0;
    for i = 1: num_phase
        ele_cell{i}(:,1) = (1:size(ele_cell{i},1))' + count;
        count = count + size(ele_cell{i},1);
    end

    % print element
    % CHEXA*,5,1,40,46,*
    % *,47
    for i = 1: num_phase
        fprintf( fid, ...
            ['CTRIA3*,%d,%d,%d,%d', ',*\n', ...
             '*,', '%d\n'], ...
            [ ele_cell{i}(:,1), i * ones(size(ele_cell{i},1),1), ele_cell{i}(:,2:4) ]' );
    end
    
    % ------------------------------------------------------------------------
	fprintf( fid, 'ENDDATA' );
    
    % ------------------------------------------------------------------------
    fclose(fid);
	
	disp('printBdf Done! Check the bdf file!');
    
end
