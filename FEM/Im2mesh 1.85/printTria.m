function printTria( vert, tria, tnum, precision_nodecoor )
% printTria: print nodes and elements into file 'test.node' & 'test.ele'.
%            Only support triangular element with 3 nodes.
%            Precision is number of digits behind decimal point, for node
%            coordinates
%            
% Revision history:
%   Jiexian Ma, mjx0799@gmail.com, Oct 2020
    
    num_node = size( vert, 1 );
    num_ele = size( tria, 1 );
    num_phase = length( unique( tnum ) );

    nodecoor_list = zeros( num_node, 3 );
    nodecoor_list( :, 1 ) = 1: num_node;
    nodecoor_list( :, 2:3 ) = vert;

    ele_list = zeros( num_ele, 5 );
    ele_list( :, 1 ) = 1: num_ele;
    ele_list( :, 2:4 ) = tria;
    ele_list( :, 5 ) = tnum;
    
    % ------------------------------------------------------------------------
    % format of number
    % field width of node numbering
    width_node_num = 1 + floor( log10( num_node ) );         % 18964 -> 5
    if width_node_num > 16
        error('more than 16 digits')
    end
    
    num_digits_of_int_part = 1 + floor( log10( max(max(vert)) ) );
                                                             % 182.9 -> 3
    if num_digits_of_int_part + precision_nodecoor + 1 > 16
        error('more than 16 digits')
    end
    format_node_coor = [ '%.', num2str( precision_nodecoor ), 'f' ];
                                                             % '%.(precision)f'
    
    % ------------------------------------------------------------------------
	fid=fopen('test.node','wW');
    % ------------------------------------------------------------------------
    fprintf( fid, '# Node count, 2d, no attribute, boundary marker off\n');
    
    fprintf( fid, '%d  2  0  0\n', num_node );
    
    fprintf( fid, '# Node index, node coordinates\n');
    
    % print node
    % 3  1.0 0.0
    fprintf( fid, ...
            [ '%d  ', format_node_coor, ' ', format_node_coor, '\n'], ...
            nodecoor_list' );
    % ------------------------------------------------------------------------
    fclose(fid);
    % ------------------------------------------------------------------------
    
    % ------------------------------------------------------------------------
	fid=fopen('test.ele','wW');
    % ------------------------------------------------------------------------
    fprintf( fid, '# Element count, nodes per element, number of attributes\n');
    
    fprintf( fid, '%d  3  %d\n', num_ele, num_phase );
    
    fprintf( fid, '# Element index, node index, attribute\n');
    % attribute is used to indicate phase
    % print elements
    % 12  2 3 4  2
    fprintf( fid, ...
            '%d  %d %d %d  %d\n', ...
            ele_list' );
    
    % ------------------------------------------------------------------------
    fclose(fid);
    % ------------------------------------------------------------------------
	
	disp('printTria Done! Check .node and .ele file!');
    
end
