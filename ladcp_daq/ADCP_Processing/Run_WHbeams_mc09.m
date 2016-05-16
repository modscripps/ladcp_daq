function ADP = Run_WHbeams_mc09(yday_beg, yday_end, WHtyp, Params)
% ADP = Run_WHbeams_mc09(yday_beg, yday_end, WHtyp, Params);
%   Converts ping beam data to velocities in geomagmetic frame,
%       returned in structure ADP .
%   Inputs:
%      yday_beg,yday_end = yearday range to process;
%      WHtyp = root name of matlab ADCP type (AD3dn,AD12up,AD12dn);
%      Params = structure with optional arguments, in fields:
%       'MatFld','MatIndx','WC_val','Psrc','Ztyp','ZGrid','Sal0','MagDec'
%          MatFld = name of directory that contains folders
%             for matlab data and matlab index (default = parent of pwd);
%          WC_val = lower correlation threshold for good beam data (def=30)

%Output
% 
% ADP = 
%            depth: [40x1 double]
%                u: [40x6481 double]
%                v: [40x6481 double]
%                w: [40x6481 double]
%             werr: [40x6481 double]
%              ec1: [40x6481 double]
%              ec2: [40x6481 double]
%              ec3: [40x6481 double]
%              ec4: [40x6481 double]
%             yday: [1x6481 double]
%           ens_no: [1x6481 double]
%         soundvel: [1x6481 double]
%        svel_calc: [1x6481 double]
%     depth_xducer: [1x6481 double]
%             temp: [1x6481 double]
%            pitch: [1x6481 double]
%             roll: [1x6481 double]
%          heading: [1x6481 double]


ADP = [];
% WHtyp must be one of the strings in cell array Mnms:
Mnms = {'SN3160m','SN3160a','SN11181m','SN7966m','SN10010m','SN11675m','SN8064m'};
iA = strmatch(WHtyp, Mnms, 'exact');
if isempty(iA)
    disp(['WHtyp = ' WHtyp ' is not one of the valid choices, exit.'])
    return
end
% For Alford mooring group, organized in ADCP s/n subfolders ...
CrName = 'mc09';
%ADtop = 'C:/mc09/mc09_adcp/Data'; % everything is be under this, by default
ADtop = '/Users/johnmickett/Cruises_Research/MC09/MC09_ADCP/Data'; %for Mickett mac (shortsands)
% Individual moored ADCPs, one for each deployment  
ADnms = {'SN3160','SN3160','SN11181','SN7966','SN10010','SN11675','SN8064'};
DeployNo = [1, 2, NaN, NaN, NaN, NaN, NaN]; % for multi deploys of same S/N
ADsub = ADnms; % if multi deploys, each has own subfolder:
for i=1:length(DeployNo)
    if ~isnan(DeployNo(i))
        ADsub{i} = fullfile(ADnms{i},['Deployment' num2str(DeployNo(i))]);
    end
end
%
% Now, define default matlab folder,index ...
MatFld = fullfile(ADtop, ADsub{iA}, 'data_mat'); % matlab data files folder
IndFld = MatFld; % matlab data index folder
MatRoot = Mnms{iA}; % for root names of matlab files and index
MatIndx = fullfile(IndFld, [MatRoot '_' CrName '_matfiles.mat']);
% and other defaults:
WC_val = 30;
%
% Assign values, special considerations for each ADCP WH deployment
EnsUP2DN = []; DepUP2DN = NaN; % for synch'd ADCPs (R-C L's EQ ones)
theta_o=20; %Angle of the ADCP beams from vertical
Cnvx = 1; % convex config
UpDown = 0; % down facing = 1, else up facing
% Pr/dbar at ADCP:
Psrc = 1; % 1=ADCP internal, 2xN vector=[yday; dbar] for interp1
Ztyp = 1; % Type of output depth grid: surface rel = 1; ADCP rel = 0;
ZGrid = [];
Sal0 = -1; % <0=use recorded soundspeed, >0(eg 35)= use for soundspeed calc
MagDec=14.1;
switch WHtyp
    case 'SN3160m'
        ZGrid = [4:4:160]'; % depth grid [m] to assign final velocities
    case 'SN3160a'
        ZGrid = [4:4:200]'; % Ascension Canyon (deploy#2)
        MagDec=14.1;
    case 'SN11181m'
        ZGrid = [8:8:300]'; % Long Ranger
    case 'SN7966m'
        ZGrid = [2:2:50]';
    case 'SN10010m'
        ZGrid = [2:2:50]';
    case 'SN11675m'
        ZGrid = [2:2:50]';
    case 'SN8064m'
        ZGrid = [32:4:220]'; % testing this unit (downward)
        UpDown = 1;
    otherwise
        disp(['Run_WHbeams_mc09:  WHtyp=' WHtyp ' is unkwown, exit'])
        return
end
if exist('Params', 'var') && isstruct(Params)
    vvs = {'MatFld','MatIndx','WC_val','Psrc','Ztyp','ZGrid','Sal0','MagDec'};
    for i=1:length(vvs)
        if isfield(Params, vvs{i}) && ~isempty(Params.(vvs{i}))
            eval([vvs{i} '=Params.(vvs{i});']);
        end
    end
end
%
% Now, get beam velocities for this ADCP
ADRaw = get_ADCP_Any(yday_beg, yday_end, MatIndx, MatFld, 1);


if isempty(ADRaw)
    disp(['Run_WHbeams_mc09: Cannot get ' WHtyp ' data, exit'])
    return
end
% Assign ADCP pressures:
if size(Psrc,2)==2 && size(Psrc,2)>1 % interpolation array
    id = find(ADRaw.yday>=Psrc(1,1) & ADRaw.yday<=Psrc(1,end));
    if ~isempty(id)
        ADRaw.depth_xducer(id) = interp1(Psrc(1,:), Psrc(2,:), ...
            ADRaw(yday(id)) );
    end
end

ADRaw.svel_calc = [];
if Sal0 > 0 % Compute actual soundspeeds
    s0 = ones(size(ADRaw.temp)) * Sal0;
    % comment out improper one (sal,pres in CU,MPa for MCG's)
    %ADRaw.svel_calc = sw_svel(s0/1000, ADRaw.temp, ADRaw.depth_xducer/100);
    ADRaw.svel_calc = sw_svel(s0, ADRaw.temp, ADRaw.depth_xducer);
    clear s0
end
%
% Convert beam data
if Ztyp > 0 % Full conversion: correct for pitch,roll, and soundspeed,
    %    map to surface-referenced depth grid and geomagnetic frame
    WHbeam_Process
    ADP=rmfield(ADP,'z_adcp');
    
else % stg1: no soundspeed correction, ADCP-relative depth bins
    WHbeam_Process_stg1
end

