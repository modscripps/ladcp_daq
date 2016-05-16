function Vel = Get_ADCP_fullSC_LADCP(RDIname, filetype, savetype, ens_ind_rng);
% Vel = Get_ADCP_fullSC(RDIname, filetype, savetype, ens_ind_rng);
%	Read in RDI file(s), put requested results in structure = Vel
% 
%  filetype: 0=process single file=RDIname, 1=multi LTA files
%		starting with RDIname (more later);
%	savetype: 0=just u,v,w, depth,time,position vectors,
%		1=more stuff (%good, etc), 2=everything
%   ens_ind_rng = min,max index values for ensembles in raw RDI file to
%       parse and convert to matlab data (11/04)
%  Vel = structure with useful RDI data
% Dave W, 03-Apr-2001, 16-nov-2004
% Get_ADCP_profRaw.m modified for raw ping data received from SWIMS ADCPs
if nargin<2
   filetype = 0; % single profile file (EG via serial Ensemble output)
end
if nargin<3
   savetype = 0;
end
if nargin<4
    ens_ind_rng = []; % get all ensembles
end

%filetype = -1; % for Raw: single file with multiple ensembles

Vel=[];

if filetype==0  % check that this is one complete file/profile
   fid=fopen(RDIname, 'r','ieee-le');
   if fid<0
      error(['Cannot open file = ' RDIname])
   end
   A = fread(fid,'char');
   frewind(fid);
   id = fread(fid,1,'uchar');
   src= fread(fid,1,'uchar');
   nbytes=fread(fid,1,'ushort');
   fclose(fid);
   if id~=127 | src~=127
      error([RDIname ' does not start with RDI header record']);
   end
   if length(A)-nbytes ~= 2
      error([RDIname ' is wrong size for single RDI profile']);
   end
end % of solo-file check

TMPfil='TMP0';

% Convert File(s)
if filetype < 1 % either: 0 = solo profile, or <0 = one file, multi-profs
   raw2mat_SC(RDIname,TMPfil, [], {'.mat'}, [], ens_ind_rng);
else
    error('Get_ADCP_fullSC.m only works for filetype<0')
end

% load in Matlab arrays, convert to useful units, save in structure
load(TMPfil)
%keyboard
delete([TMPfil '.mat']);
% Check if no valid data were found, exit if so
if isempty(ensemble_number)
    Vel = [];
    return
end

BADVEL = -32768;

UpDown = 1; % -1=down, +1=upward (% upward for SC in aeg04)

Vel.ens_no = ensemble_number + ensMSB*65536;

% Some will be the same for all ensembles:
Vel.z_adcp =  - UpDown*(dis1/100 + [0:nbins-1]*binlen/100); % meters
Vel.p_adcp = Vel.z_adcp/100; % MPa
Vel.pulselen = pulselen;

% Others will vary:
Vel.yday = yearday(day,month,year,hour,minute,second+hundreths/100);
% %% BT - range and velocity
% if isempty(btrange1) % If set for NOT bottom-tracking, put in NaNs (3-18-2003)
%         ix = NaN*Vel.ens_no;
%         btrange1 = ix; btrange2 = ix; btrange3 = ix; btrange4 = ix;
%         btvel1 = ix; btvel2 = ix; btvel3 = ix; btvel4 = ix;
% end

% Orientation data (Internal for SWIMS)
Vel.heading = heading/100;
Vel.pitch = pitch/100;
Vel.roll = roll/100;
Vel.hdg_std = stdhed/1;
Vel.pitch_std = stdpitch/10;
Vel.roll_std = stdroll/10;
Vel.depth_xducer = xducerdepth/10;
Vel.soundvel = soundspeedRDI;
Vel.temp = degC/100;

if exist('btrange1')
btrange = [btrange1; btrange2; btrange3; btrange4];
for i=1:4
    ib = find(btrange(i,:) == 0);
    if ~isempty(ib)
        btrange(i,ib) = NaN; 
    end
end
Vel.btrange = btrange;
Vel.bottomBT = nanmean(btrange)/100 + xducerdepth/10;
end

if exist('btvel1')
btvel1=btvel1;btvel2=btvel2;btvel3=btvel3;btvel4=btvel4;
ib = find(btvel1==BADVEL);
if ~isempty(ib)
    btvel1(ib) = NaN; 
end
ib = find(btvel2==BADVEL);
if ~isempty(ib)
    btvel2(ib) = NaN;
end
ib = find(btvel3==BADVEL);
if ~isempty(ib)
    btvel3(ib) = NaN;
end
Vel.btvel_bm = [btvel1; btvel2; btvel3; btvel4] / 1000; % m/s
end

%% Measured beam components
bad = find(v1==BADVEL); if bad, v1(bad) = NaN; end
bad = find(v2==BADVEL); if bad, v2(bad) = NaN; end
bad = find(v3==BADVEL); if bad, v3(bad) = NaN; end
bad = find(v4==BADVEL); if bad, v4(bad) = NaN; end
% These are in earth's frame for typical SC deployment:
Vel.u_wat = v1/1000; Vel.v_wat = v2/1000; % m/s
Vel.w_wat = v3/1000; Vel.err_wat = v4/1000; % m/s
        
if savetype>0 
    % Save more ...
    Vel.ec1_bm = e1; Vel.ec2_bm = e2; Vel.ec3_bm = e3; Vel.ec4_bm = e4;
    Vel.cor1_bm = cor1; Vel.cor2_bm = cor2; Vel.cor3_bm = cor3; Vel.cor4_bm = cor4;
    Vel.pg1 = pg1; Vel.pg2 = pg2; Vel.pg3 = pg3; Vel.pg4 = pg4;
end

return



