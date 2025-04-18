function printInp_multiPart(name, nodecoor_cell, ele_cell, ele_type, precision_nodecoor )
%printInp: print the nodes and elements into Inp file 'test.inp',
%          test in software Abaqus. Each phase corresponds to one part in
%          Abaqus.
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
	% Parts
	%   
	% Assembly
	%	Instance
    %
    % 	Nset
    % 		xmin, xmax, ymin, ymax
    % 		interfnode
	%	Surface
    % 	
	%	Constraint
    % ------------------------------------------------------------------------
	
    
    % get node set, used to define interface and boundary condition
    interfnode_cell = getInterf( nodecoor_cell );
    [ xmin_node_cell, xmax_node_cell, ...
      ymin_node_cell, ymax_node_cell ] = getBCNode( nodecoor_cell );
  
    % ------------------------------------------------------------------------
    % initialize
    num_phase = length( nodecoor_cell );
    part_ascii = 65: ( 65 + num_phase - 1);
    part_char = char( part_ascii );     % 'ABCD...'
    
    % node numbering of last node in each nodecoor_cell{i}
    end_node_v = cellfun( @(A) A(end,1), nodecoor_cell ); 
    num_node = max(end_node_v);
    
    % number of element in each ele_cell{i}
    [ele_size_v, ~] = cellfun( @size, ele_cell );
    num_ele = sum(ele_size_v);
    
    % ------------------------------------------------------------------------
    % format of number
    % field width of node numbering
    width_node_num = 1 + floor( log10( num_node ) );         % 18964 -> 5
    % print format of node numbering
    format_node_num = [ '%', num2str(width_node_num), 'd' ]; % '%5d'
    
    max_coor = max( cellfun( @(A) max(max(A(:,2:3))), nodecoor_cell ) );
    num_digits_of_int_part = 1 + floor( log10( max_coor ) );
                                                             % 182.9 -> 3
    % print format of node coordinates
    format_node_coor = [ '%', ...
                num2str( num_digits_of_int_part + precision_nodecoor + 1 ), '.', ...
                num2str( precision_nodecoor ), 'f' ];
                                          % '%(3+1+precision).(precision)f'
    
    % field width of element numbering
    width_ele_num = 1 + floor( log10( num_ele ) );          % 18964 -> 5
    % print format of element numbering
    format_ele_num = [ '%', num2str(width_ele_num), 'd' ];  % '%5d'
    
    % ------------------------------------------------------------------------
	fid=fopen([name, '.inp'],'wW');
    % ------------------------------------------------------------------------
	% Heading
    fprintf( fid, [...
        '*Heading'                                              '\n'...
        '*Preprint, echo=NO, model=NO, history=NO, contact=NO'  '\n'...
        ] );
	% ------------------------------------------------------------------------
	% Parts
    for i = 1: num_phase
        fprintf( fid, [...
            '** PARTS'                  '\n'...
            '*Part, name=Part-%c'       '\n'...
            '*Node'                     '\n'...
            ], part_char(i) );
        
        fprintf( fid, ...
            [ format_node_num, ',', format_node_coor, ',', ...
                                    format_node_coor, '\n'], ...
            nodecoor_cell{i}' );
        
        fprintf( fid, [...
            '*Element, type=%s'  '\n'...
            ], ele_type );
        printEle( fid, ele_cell{i}, format_ele_num, format_node_num );
        
        fprintf( fid, '*Nset, nset=Set-%c\n', part_char(i) );
        printSet( fid, nodecoor_cell{i}(:,1) );
        
        fprintf( fid, '*Elset, elset=Set-%c\n', part_char(i) );
        printSet( fid, ele_cell{i}(:,1) );
        
        fprintf( fid, [...
            '** Section: Section-%c'            '\n'...
            '*Solid Section, elset=Set-%c, material=Material-%c'  '\n'...
            ','                                 '\n'...
            '*End Part'                         '\n'...
            '**'                                '\n'...
            ], ...
            part_char(i), ...
            part_char(i), part_char(i) );
    end
	
	% ------------------------------------------------------------------------
	% Assembly
    fprintf( fid, [...
        '** ASSEMBLY'               '\n'...
        '*Assembly, name=Assembly'	'\n'...
        '**'                        '\n'...
        ] );
	% ------------------------------------------------------------------------
	
    % Instance ABCD...
    for i = 1: num_phase
        fprintf( fid, [...
            '*Instance, name=Part-%c-1, part=Part-%c'   '\n'...
            '*End Instance'                     '\n'...
            '**'                                '\n'...
            ], ...
            part_char(i), part_char(i) );
    end
    
    % ------------------------------------------------------------------------
    % Nset
    
    % xmin
    for i = 1: num_phase
        if ~isempty( xmin_node_cell{i} )
            fprintf( fid, [...
                '*Nset, nset=Set-%c-Xmin, instance=Part-%c-1'   '\n'...
                ], part_char(i), part_char(i) );
            printSet( fid, xmin_node_cell{i} );
        end
    end
    
    % xmax
    for i = 1: num_phase
        if ~isempty( xmax_node_cell{i} )
            fprintf( fid, [...
                '*Nset, nset=Set-%c-Xmax, instance=Part-%c-1'   '\n'...
                ], part_char(i), part_char(i) );
            printSet( fid, xmax_node_cell{i} );
        end
    end
    
    % ymin
    % u can use this node set to define boundary condition
    % example:
    % ** Name: BC-FixY Type: Displacement/Rotation
    % *Boundary
    % Set-A-Ymin, 2, 2
    for i = 1: num_phase
        if ~isempty( ymin_node_cell{i} )
            fprintf( fid, [...
                '*Nset, nset=Set-%c-Ymin, instance=Part-%c-1'   '\n'...
                ], part_char(i), part_char(i) );
            printSet( fid, ymin_node_cell{i} );
        end
    end
    
    % ymax
    for i = 1: num_phase
        if ~isempty( ymax_node_cell{i} )
            fprintf( fid, [...
                '*Nset, nset=Set-%c-Ymax, instance=Part-%c-1'   '\n'...
                ], part_char(i), part_char(i) );
            printSet( fid, ymax_node_cell{i} );
        end
    end
    
    % node at interface
    for i = 1: num_phase-1
        for j = i+1: num_phase
            if ~isempty( interfnode_cell{i,j} )
                % node in part i of interface i,j
                fprintf( fid, [...
                    '*Nset, nset=Set-%c-Interf%c%c, instance=Part-%c-1' '\n'...
                    ], part_char(i), part_char(i), part_char(j), part_char(i) );
                printSet( fid, interfnode_cell{i,j} );
                
                % node in part j of interface i,j
                fprintf( fid, [...
                    '*Nset, nset=Set-%c-Interf%c%c, instance=Part-%c-1' '\n'...
                    ], part_char(j), part_char(i), part_char(j), part_char(j) );
                printSet( fid, interfnode_cell{j,i} );
            end
        end
    end
    
	fprintf( fid, '%s\n', '**' );

    % ------------------------------------------------------------------------
    % surface (define surface using node sets)
    % Note: node based surface is not supported in abaqus CAE
    
    % ymax surface
    % u can use Inp file to define the Interaction between ymax and punch
    % example:
    % ** Interaction: Int-PunchSpecimen
    % *Contact Pair, interaction=IntProp-NoFrition, type=SURFACE TO SURFACE, tracking=STATE
    % Set-B-Ymax_CNS_, Punch-1.Surf-Punch
    
    for i = 1: num_phase
        if ~isempty( ymax_node_cell{i} )
            fprintf( fid, [...
                '*Surface, type=NODE, name=Set-%c-Ymax_CNS_, internal'  '\n'...
                'Set-%c-Ymax, 1.'     '\n'...
                ], ...
                part_char(i), part_char(i) );
        end
    end
     
    % interface
    for i = 1: num_phase-1
        for j = i+1: num_phase
            if ~isempty( interfnode_cell{i,j} )
                % node in part i of interface i,j
                fprintf( fid, [...
                    '*Surface, type=NODE, name=Set-%c-Interf%c%c_CNS_, internal' '\n'...
                    'Set-%c-Interf%c%c, 1.'     '\n'...
                    ], ...
                    part_char(i), part_char(i), part_char(j), ...
                    part_char(i), part_char(i), part_char(j) );
                
                % node in part j of interface i,j
                fprintf( fid, [...
                    '*Surface, type=NODE, name=Set-%c-Interf%c%c_CNS_, internal' '\n'...
                    'Set-%c-Interf%c%c, 1.'     '\n'...
                    ], ...
                    part_char(j), part_char(i), part_char(j), ...
                    part_char(j), part_char(i), part_char(j) );                
            end
        end
    end
    
    % ------------------------------------------------------------------------
	% Constraint
	fprintf( fid, '%s\n', '** Constraint' );
    
    for i = 1: num_phase-1
        for j = i+1: num_phase
            if ~isempty( interfnode_cell{i,j} )
                fprintf( fid, [...
                    '*Tie, name=Constraint-Interf%c%c, adjust=yes, type=NODE TO SURFACE' '\n'...
                    'Set-%c-Interf%c%c_CNS_, Set-%c-Interf%c%c_CNS_'  '\n'...
                    ], ...
                    part_char(i), part_char(j), ...
                    part_char(i), part_char(i), part_char(j), ...
                    part_char(j), part_char(i), part_char(j) );
            end
        end
    end
    
	% ------------------------------------------------------------------------
	fprintf( fid, '%s', '*End Assembly' );
	
    fclose(fid);
	
	disp('printInp_multiPart Done! Check the inp file!');
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

function printSet( fid, nodeset )

    for i=1:length(nodeset)
        if mod( i, 16 ) == 0 || i == length(nodeset)
            fprintf( fid, '%d\n', nodeset(i) );
        else
            fprintf( fid, '%d, ', nodeset(i) );
        end
    end
end