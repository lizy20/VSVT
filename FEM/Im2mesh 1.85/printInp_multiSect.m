function printInp_multiSect( nodecoor_list, ele_cell, ele_type, precision_nodecoor )
% printInp_multiSect: print the nodes and elements into Inp file 'test.inp',
%          test in software Abaqus. One part with multiple sections. Each
%          phase corresponds to one section in Abaqus.
%
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, Oct 2020
% Cite As
%   Jiexian Ma (2020). Im2mesh (2D image to triangular meshes) (https://ww
%   w.mathworks.com/matlabcentral/fileexchange/71772-im2mesh-2d-image-to-t
%   riangular-meshes), MATLAB Central File Exchange. Retrieved May 21, 202
%   0.
%   

	% format of Inp file
	% ------------------------------------------------------------------------
	% Heading
	%
	% Node
    %
    % Element
    %
    % Section
    % ------------------------------------------------------------------------
	
    num_sect = length( ele_cell );
    sect_ascii = 65: ( 65 + num_sect - 1);
    sect_char = char( sect_ascii );     % 'ABCD...'
    
    % format of number
    % '%.(precision)f'
    format_node_coor = [ '%.', num2str( precision_nodecoor ), 'f' ];
    
    format_node_num = '%d';
    format_ele_num = '%d';
    
    % ------------------------------------------------------------------------
	fid=fopen('test_sections.inp','wW');
    % ------------------------------------------------------------------------
	% Heading
    fprintf( fid, [...
        '*Heading'                                              '\n'...
        '*Preprint, echo=NO, model=NO, history=NO, contact=NO'  '\n'...
        '**'                                                    '\n'...
        ] );
    
	% ------------------------------------------------------------------------
    % Node
    fprintf( fid, '*Node\n' );
    
    % print coordinates of nodes
    % '%d,%.4f,%.4f,%.4f\n'
    fprintf( fid, ...
        [ format_node_num, ',', format_node_coor, ',', ...
                                format_node_coor, '\n'], ...
        nodecoor_list' );
    
    % ------------------------------------------------------------------------
    % Element
    for i = 1: num_sect
        fprintf( fid, [...
            '*Element, type=%s, elset=Set-%c'  '\n'...
            ], ele_type, sect_char(i) );
        
        printEle( fid, ele_cell{i}, format_ele_num, format_node_num );
    end
    
    % ------------------------------------------------------------------------
    % Section
    for i = 1: num_sect
        fprintf( fid, [...
            '** Section: Section-%c'            '\n'...
            '*Solid Section, elset=Set-%c, material=Material-%c'  '\n'...
            ','                                 '\n'...
            ], ...
            sect_char(i), ...
            sect_char(i), sect_char(i) );
    end
    
	fprintf( fid, '**' );
    
    % ------------------------------------------------------------------------
    fclose(fid);
	
	disp('printInp_multiSect Done! Check the inp file!');
end


function printEle( fid, ele, format_ele_num, format_node_num )
% work for linear element and quadratic element

    num_column = size( ele, 2 );
    switch num_column
    case 4
        % linear element
        fprintf( fid, ...
            [ format_ele_num, ',', repmat([format_node_num, ','], [1,2]), format_node_num, '\n' ], ...
            ele' );
        
    case 7
        % quadratic element
        fprintf( fid, ...
            [ format_ele_num, ',', repmat([format_node_num, ','], [1,5]), format_node_num, '\n' ], ...
            ele' );
        
    otherwise
        disp('unidentified data')
    end

end