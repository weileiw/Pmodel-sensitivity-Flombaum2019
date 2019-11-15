load omega2_ad1e-05_ai1000.mat % load transport data
load po4obs_91x180x24.mat      % WOA2013 DIP 
load npp_91x180.mat npp        % Satellite NPP

ineg = find(po4obs(:)<0); % find and reset negative DIP;
po4obs(ineg) = 0.02;

TRdiv = -output.TR;         % transport divergence [s^-1];
M3d   = output.M3d;         % 3D land-sea mask;
grd   = output.grid;        % grid data
dVt   = grd.DXT3d.*grd.DYT3d.*grd.DZT3d; % 3D grid volume;
% prescribe some constants.
dpa = 365;                  % days per year;
spd = 24*60^2;              % seconds per day;
spa = dpa*spd;              % seconds per year;
iwet = find(M3d(:)==1);     % wet point index;
nwet = length(iwet);
I = speye(nwet);

parm.p2c = (1/106)*M3d;     % carbon to phosphate ratio

parm.sigma = 1/3;             % portion of production to DOP[unitless]
parm.kappa_p = 1/(720*60^2);   % POP disolution constant [s^-1];
parm.kappa_g = 1/(1e6*spa);    % DIP geological restoring constant [s^-1];

% global averaged DIP conc. [mmol m^-3];
parm.DIPbar  = sum(po4obs(iwet).*dVt(iwet))./sum(dVt(iwet));

inan = find(isnan(npp(:))); npp(inan) = 0; % get rid of invalid npp;
parm.Lambda = M3d*0;

parm.Lambda(:,:,1) = 0.5*(1/grd.dzt(1))*parm.p2c(:,:,1)./(1e-9+po4obs(:,:,1));
parm.Lambda(:,:,2) = 0.5*(1/grd.dzt(2))*parm.p2c(:,:,2)./(1e-9+po4obs(:,:,2));
parm.Lambda(:,:,3:end) = 0;
npp = npp/(12*spd);
parm.npp = npp;

load xhat_91x180_control.mat % load optimimal parameters;
p = R.xhat;

% coefficients to adjust b and kappa_d
%co = [0.85,0.90,0.95,1.05,1.10,1.15]; 
% sigma = [0.10,0.20,0.40,0.50,0.60,0.70];
sigma = 0.2;
for ji = 1:length(sigma)
    
    % parm.b = p(1)*co(ji);
    % parm.kappa_d = p(2)*co(ji);
    parm.sigma = sigma(ji);
    P = eqPcycle_sensi(parm,grd,M3d,TRdiv,p);
    
    DIP = M3d+nan;
    POP = M3d+nan;
    DOP = M3d+nan;
    
    DIP(iwet) = P(1:nwet);
    POP(iwet) = P(nwet+1:2*nwet);
    DOP(iwet) = P(2*nwet+1:end);
    
    fname = sprintf('R91x180_whole_sigma_%2d_percent',sigma(ji)*100);
    save(fname,'DIP','POP','DOP');
    
end
