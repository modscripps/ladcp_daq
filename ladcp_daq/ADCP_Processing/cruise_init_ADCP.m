% cruise_init_ADCP.m  -  initialize folders, indices for cruise;
%	Edit values above line='%%%%%%%%%%%%%%%%%% HERE %%%%%%%%%%%%%%%%%%'
%   for new cruise.  cd to X:/swims/'cruisename' before running this.
% Dave W, apr-2002
% THIS one just for doing bb150 files for may 2006 tests

CR_name = 'IWISE10';
SWIM_no = 0;
% structures for matlab index files:
Cruise.name = 'IWISE10';
Cruise.start_date = '14-Aug-2010 00:00:00';  % UTC, approx to start
Cruise.end_date = '11-Sep-2010 00:00:00';
[ydayS,yr]=datenum2yday(datenum(Cruise.start_date));
ydayE,yr]=datenum2yday(datenum(Cruise.start_date));

Def(1).year = yr;
Def(1).start_yday = ydayS;
Def(1).end_yday = ydayE;
Def(1).yday_seconds_offset = 0;
Def(1).SWIMS_num = SWIM_no;
Idx(1).yday_beg = [];
Idx(1).yday_end = [];
Idx(1).filename = [];
PROG ='cruise_init_ADCP.m';

cr_dir=pwd;
DATA_dir = '/Users/johnmickett/Cruises_Research/IWISE10/ADCP/Data/data_mat';

%%%%%%%%%%%%%%%%%% HERE %%%%%%%%%%%%%%%%%%

data_typs={'ADCP'};

% then matlab (stage 1) index; 
for id=1:length(data_typs)
    if id<=length(data_typs)
        fnidx = [data_typs{id} '_' CR_name '_matfiles.mat'];
        disp(['Init ' fnidx])
        if exist(fullfile('.',fnidx), 'file') & 0
            disp('   Already exists, skipping !!')
        else
            clear Set_params Index
            Set_params = Def;  Index = Idx;
            % add more fields to Set_params(1), Index(1) if needed
            %
            save( fullfile('.',fnidx), ...
                'Cruise','Set_params','Index','PROG');
        end
    end % of initializing index into Matlab (stage 1) files
end