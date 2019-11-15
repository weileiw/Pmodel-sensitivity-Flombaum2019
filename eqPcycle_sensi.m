function [P] = eqPcycle_sensi(parm,grd,M3d,TRdiv,p)
% output: P is model prediction of DIP,POP,and DOP
% output: D is dPdp, which is used for parameter optimization
% output: D2 is d2Pdp2, which is used for parameter optimization
% M*P = r; dMdp*P+M*dPdp = drdp; dPdp = M\(drdp-dMdp*P);
    
    iocn   = find(M3d(:));        % wet point index;
    n_iocn = length(iocn);        % number of wet points;
    I = speye(n_iocn);            % make an identity matrix;
    % fixed parameters
    DIPbar  = M3d(iocn)*parm.DIPbar;  % gobal arerage PO4 conc.[mmol m^-3]; 
    kappa_g = parm.kappa_g;           % PO4 geological restore const.[s^-1];
    kappa_p = parm.kappa_p;           % POP solubilization rate constant
    sigma   = parm.sigma;             % portion of production to DOP[unitless]
    
    
    b = p(1); % Martin curve exponent;
    kappa_d = p(2);  % DOP remineralization const.[s^-1];
    alpha   = p(3);  % npp scaling factor for DIP uptake rate
    beta    = p(4);  % npp scaling exponent for DIP uptake rate
    
    PFdiv = buildPFD(M3d,b,kappa_p,grd); % particle flux divergence operator;
    
    npp     = parm.npp;               % net primary production 
    Lambda  = parm.Lambda;    
    Lambda(:,:,1) = (npp.^beta).*Lambda(:,:,1);
    Lambda(:,:,2) = (npp.^beta).*Lambda(:,:,2);
    L = d0(Lambda(iocn));  % PO4 assimilation rate [s^-1];
    
    % build Jacobian matrix
    % +++++++++++++++++++++++++++++++++++++++++
    % column 1 (d/dDIP)
    J{1,1} = TRdiv+alpha*L+kappa_g*I;
    J{2,1} = -(1-sigma)*alpha*L;
    J{3,1} = -sigma*alpha*L;
    
    % column 2 (d/dPOP)
    J{1,2} = 0*I;
    J{2,2} = PFdiv+kappa_p*I;
    J{3,2} = -kappa_p*I;
    
    % column 3 (d/dDOP)
    J{1,3} = -kappa_d*I;
    J{2,3} = 0*I;
    J{3,3} = TRdiv+kappa_d*I;
    % ++++++++++++++++++++++++++++++++++++++++
    
    % right hand side of phosphate equations 
    RHS{1,1} = DIPbar*kappa_g;
    RHS{2,1} = sparse(n_iocn,1);
    RHS{3,1} = sparse(n_iocn,1);
    
    % factoring Jacobian matrix
    Jac = cell2mat(J);
    FJ = mfactor(Jac); 
    % solve solutions 
    P = mfactor(FJ,cell2mat(RHS));
    
end

% --------------------------------------------------------------------

function [PFdiv,dPFDdb,d2PFDdb2,dPFDdkappa_p] ...
        = buildPFD(M3d,b,kappa_p,grd);
    
    [ny,nx,nz] = size(M3d);
    M3D = zeros(ny,nx,nz+1);
    M3D(:,:,1:end-1) = M3d;
    % add the zw coordinate at the top of the extra layer
    ZW3d = grd.ZW3d;
    ZW3d = ZW3d(:,:,[1:end,end]);
    ZW3d(:,:,end) = grd.ZW3d(:,:,end)+grd.dzt(end);
    % areas of the top of the grid box
    dAt = (grd.DXT3d.*grd.DYT3d).*M3d;
    % volume of the grid boxes
    dVt = (grd.DXT3d.*grd.DYT3d.*grd.DZT3d).*M3d;
    %
    n = nx*ny*(nz+1);
    I0 = speye(n);
    i0 = zeros(ny,nx,nz+1);
    i0(:) = 1:n;
    % periodic shifts OK because M3D has a layer of zeros on the bottom
    iu = i0(:,:,[nz+1,1:nz]); %use a periodic upward shift
    ib = i0(:,:,[2:nz+1,1]); % use a periodic downward shift
    IU = I0(iu,:);
    IB = I0(ib,:);
    % keep only wet boxes
    iocn = find(M3D(:));
    I0 = I0(iocn,:); I0 = I0(:,iocn);
    IU = IU(:,iocn); IU = IU(iocn,:);
    IB = IB(:,iocn); IB = IB(iocn,:);
    % (averages POP onto the top of the grid boxes)
    AVG = d0((I0+IU)*M3D(iocn))\(I0+IU);
    % (compute the divergence in the center of the grid boxes)
    DIV = d0(dVt(iocn))\(I0-IB)*d0(dAt(iocn));
    % (compute the flux at the top of the grid cells)
    % mimics a Martin curve flux attenuation profile
    %(see Kriest and Oschelies 2008 in Biogeosciences)
    r = kappa_p;
    a = r/b;
    % particle sinking velocity at the top of the grid cells.
    MSK = M3D.*M3D(:,:,[nz+1,1:nz]);
    M = MSK.*ZW3d;
    w = -a*M;
    dadb = -r/(b^2);
    dwdb = -dadb*M;
    d2adb2 = 2*r/(b^3);
    d2wdb2 = -d2adb2*M;
    dadr = 1./b;
    dwdr = -dadr*M;
    %FLUX = d0(w(iocn))*AVG;
    FLUX = d0(w(iocn))*IU;
    dFLUXdb = d0(dwdb(iocn))*IU;
    d2FLUXdb2 = d0(d2wdb2(iocn))*IU;
    dFLUXdr = d0(dwdr(iocn))*IU;
    % particle flux divergence operator
    PFdiv = DIV*FLUX;
    dPFDdb = DIV*dFLUXdb;
    d2PFDdb2 = DIV*d2FLUXdb2;
    dPFDdkappa_p = DIV*dFLUXdr;    
end