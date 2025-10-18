function StepwiseCMD_demo
% Demonstration for multivariate 3-D random field simulation using stepwise 
%   covariance matrix decomposition method
% Written by Te XIAO [xiaote@sjtu.edu.cn] (2018/12/01)
% For more details, please refer to:
%     Li DQ, Xiao T, Zhang LM, & Cao ZJ. (2019). Stepwise covariance 
%     matrix decomposition for efficient simulation of multivariate 
%     large-scale three-dimensional random fields. Applied Mathematical 
%     Modelling, 68, 169-181. DOI: 10.1016/j.apm.2018.11.011

clc; clear; close all;

% Scale of fluctuation
tx = 30; ty = 20; tz = 1;
% Separable correlation function
CorrFun1 = @(d, t) exp(-2*d/t);
% Node domain, dm = [min, max]
dmx = [0, 100]; dx = 1;
dmy = [0, 100]; dy = 1;
dmz = [0, 50];  dz = 1;
% Node
x = dmx(1):dx:dmx(2); nx = length(x);
y = dmy(1):dy:dmy(2); ny = length(y);
z = dmz(1):dz:dmz(2); nz = length(z);
% Relative distance
[x1, x2] = ndgrid(x); Dx = abs(x1-x2);
[y1, y2] = ndgrid(y); Dy = abs(y1-y2);
[z1, z2] = ndgrid(z); Dz = abs(z1-z2);
% Correlation matrix & Cholesky decomposition
Rx = CorrFun1(Dx, tx); Lx = chol(Rx, 'lower');
Ry = CorrFun1(Dy, ty); Ly = chol(Ry, 'lower');
Rz = CorrFun1(Dz, tz); Lz = chol(Rz, 'lower');
% Cross-correlation matrix
np = 3;
Rc = eye(np);
Rc(1, 2) = 0.8; Rc(2, 1) = Rc(1, 2);
Rc(1, 3) = -0.8; Rc(3, 1) = Rc(1, 3);
Rc(2, 3) = -0.8; Rc(3, 2) = Rc(2, 3);
Lc = chol(Rc, 'lower');

% Independent random variable
U = randn(nx, ny, nz, np);        % nx-ny-nz-np array
% Random field simulation, U --> X
X = Lx*reshape(U, nx, ny*nz*np);  % nx-ny*nz*np matrix
X = reshape(X, nx, ny, nz, np);   % nx-ny-nz-np array
X = permute(X, [2, 3, 4, 1]);     % ny-nz-np-nx array
X = Ly*reshape(X, ny, nz*np*nx);  % ny-nz*np*nx matrix
X = reshape(X, ny, nz, np, nx);   % ny-nz-np-nx array
X = permute(X, [2, 3, 4, 1]);     % nz-np-nx-ny array
X = Lz*reshape(X, nz, np*nx*ny);  % nz-np*nx*ny matrix
X = reshape(X, nz, np, nx, ny);   % nz-np-nx-ny array
X = permute(X, [2, 3, 4, 1]);     % np-nx-ny-nz array
X = Lc*reshape(X, np, nx*ny*nz);  % np-nx*ny*nz matrix
X = reshape(X, np, nx, ny, nz);   % np-nx-ny-nz array
X = permute(X, [2, 3, 4, 1]);     % nx-ny-nz-np array

% Random field visualization
cmin = -4; cmax = 4;
X = reshape(X, nx*ny*nz, np);     % nx*ny*nz-np matrix
[xt, yt, zt] = meshgrid(x, y, z); % ny-nx-nz array
for ip = 1:np
    Xi = reshape(X(:, ip), nx, ny, nz);
    Xi = permute(Xi, [2, 1, 3]); % ny-nx-nz array
    figure; colormap('jet');
    hs = slice(xt, yt, zt, Xi, dmx, dmy, dmz);
    set(hs, 'EdgeAlpha', 0.1);
    caxis([cmin, cmax]); axis equal;
    xlabel('x'); ylabel('y'); zlabel('z');
end
return