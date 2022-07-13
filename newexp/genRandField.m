%%%
%%% genRandField.m
%%%
%%% Generates a random two-dimensional field. Useful for creating random 
%%% initial conditions or random topography. The output field 'psi' will
%%% have characteristic spectral wavelength 'lambda' with an exponential 
%%% width 'W' in spectral space, and and rms amplitude 'psirms'.
%%% Nx and Ny define the grid size, and Lx and Ly define the 
%%% zonal and meridional domain lengths. If 'W' is an empty vector then a
%%% default value of 1/8 * the wavenumber corresponding to lambda will be
%%% used.
%%% 
function F = genRandField (lambda,W,Frms,Nx,Ny,Lx,Ly) 
 
  %%% Spectral grids  
  k = [0:1:Nx/2-1,-Nx/2:1:-1]; %%% Zonal wavenumber
  K_xk = 2*pi.*(k)./Lx;
  l = [0:1:Ny/2-1,-Ny/2:1:-1]; %%% Meridional wavenumber
  K_yl = 2*pi.*(l)./Ly;
  [K_ykl,K_xkl] = meshgrid(K_yl, K_xk);   
  
  %%% Most energetic wavenumber
  K_0 = 2*pi/lambda; 
  
  %%% Exponential width of energy band in wavenumber space
  if (isempty(W))
    W = K_0/8; 
  end

  %%% Amplitude is exponential about K0, and phases are random. N.B. here
  %%% we only define the amplitude up to a constant - below we constrain it.  
  K = sqrt(K_xkl.^2 + K_ykl.^2);
  theta = 2 .* pi .* rand(Nx,Ny);
  psi_fft = K.^(-1).*exp(-((K-K_0)/W).^2) .* exp(1i*theta);

  %%% Avoids infinite mode-0 amplitude 
  psi_fft(1,1) = 0;

  %%% Transform back to real space
  F = real(ifft2(psi_fft));

  %%% Normalize so that RMS is Frms
  F = F * Frms./sqrt(sum(sum(F.^2))/(Nx*Ny));  

end
