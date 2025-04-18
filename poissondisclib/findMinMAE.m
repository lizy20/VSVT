
function out = findMinMAE( est, truth, varargin )
  % Result is min_k MAE( k * est - truth )
  %
  % out = findMinMAE( est, truth [, 'verbose', true/false ] )
  %
  % Written by Nicholas Dwork, Copyright 2024
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 2
    disp( 'Usage:  out = findMinMAE( recon, trueRecon [, ''verbose'', true/false ] )' );
    if nargout > 0, out = []; end
    return;
  end

  p = inputParser;
  p.addParameter( 'verbose', false );
  p.parse( varargin{:} );
  verbose = p.Results.verbose;

  f = @(k) mae( k * est, truth );

  LB = 0;
  UB = max( est(:) ) / mean( truth(:) ) * 10;

  k = goldenSectionSearch( f, LB, UB, 'tol', 1d-6, 'verbose', verbose );
  out = f( k );
end
