
function bp = ctBackProject( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
  dx, dy, varargin )
  % bp = ctBackProject( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy [, 'type', ] )
  %
  % Performs back projection of a sinogram with physical units of a computed
  % tomography detector
  %
  % Inputs:
  % sinogram - 2D array, each row corresponds to a theta, each column to a detector
  % thetas - a 1D array of thetas corresponding to each sinogram row
  % dSize - the size of each detector in meters
  % cx - center x position of reconstructed region
  % cy - center y position of reconstructed region
  % xFOV - horizontal field of view in meters
  % yFOV - vertical field of view in meters
  % Nx - number of pixels horizontally in reconstruction
  % Ny - number of pixels vertically in reconstruction
  % dx - horizontal size of pixel in meters
  % dy - vertical size of pixel in meters
  % type (optional) - 'iso' (default) or 'fast'
  %    'iso' uses a rotation that's an isometry
  %    'fast' faster implementation
  %
  % Written by Nicholas Dwork - Copyright 2013
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  defaultType = 'fast';
  expectedTypes = { 'iso', 'fast' };

  p = inputParser;
  p.addRequired( 'sinogram', @(x) ismatrix(x) );
  p.addRequired( 'thetas', @(x) ismatrix(x) );
  p.addRequired( 'dSize', @isnumeric );
  p.addRequired( 'cx', @isnumeric );
  p.addRequired( 'cy', @isnumeric );
  p.addRequired( 'Nx', @isnumeric );
  p.addRequired( 'Ny', @isnumeric );
  p.addRequired( 'dx', @isnumeric );
  p.addRequired( 'dy', @isnumeric );
  p.addOptional( 'type', defaultType, @(x) any( validatestring(x,expectedTypes) ) );

  if dx ~= dy
    error( 'ctBackProject: dx and dy must be equal for this version of back projection' );
  end

  [nThetas, nDetectors] = size(sinogram);

  dOffset = 0;   % detector center offset
  dLocs = ( (0:nDetectors-1) - floor(0.5*nDetectors) ) * dSize - dOffset;

  % create the adjoint of the interpolation matrix
  if mod( Nx, 2 )==0
    locs = ( ( 0 : Nx-1 ) - 0.5 * Nx + 0.5 ) * dx;
  else
    locs = ( ( 0 : Nx-1 ) - floor( 0.5 * Nx ) ) * dx;
  end
  B = zeros(nDetectors, Nx);
  for i = 1 : Nx
    tmp = zeros(1,Nx);
    tmp(i) = 1;
    interped = interp1( locs, tmp, dLocs, 'linear', 0 );
    B(:,i) = interped;
  end
  Bt = B';

  bp = zeros( Ny, Nx );
  parfor thIndx = 1 : nThetas
    theta = thetas( thIndx );
    line = sinogram( thIndx, : );
    tmp = Bt * line';
    smeared = ( dy * ones(Ny,1) ) * tmp';
    bp = bp + rotImg( smeared, theta, 'op', 'transp' );
    if mod( thIndx, 5 ) == 0
      disp([ 'ctBackProject Theta Indx / Theta: ', num2str(thIndx), '/', num2str(theta) ]);
    end
  end

end
