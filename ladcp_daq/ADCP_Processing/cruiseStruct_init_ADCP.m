% cruiseStruct_init_ADCP.m  -  initialize folders, indices for cruise;
%cd to the mfiles directory /cruisename/ADCP/mfiles before running this!
% Dave W, apr-2002, JBM Jul 2009


%this is called by Conv_ADCP_mooring_CRUISENAME.m

CR_name = InP.Cruise.Name;
%SWIM_no = 0;
% structures for matlab index files:
Cruise.name = InP.Cruise.Name;
Cruise.start_date = InP.Cruise.StartTime;  % UTC, approx to start
Cruise.end_date = InP.Cruise.EndTime;

%VV=datevec(InP.Cruise.StartTime);
[ydayy,yrr]=datenum2yday(datenum(InP.Cruise.StartTime));
%[ydayy2,yrr2]=datenum2yday(datenum(InP.Cruise.EndTime));


ydayy2=ydayy+[datenum(InP.Cruise.EndTime)-datenum(InP.Cruise.StartTime)];



Def(1).start_year = yrr;
Def(1).start_dtnum = datenum(InP.Cruise.StartTime); 
Def(1).end_dtnum = datenum(InP.Cruise.EndTime); %or 
Def(1).dtnum_seconds_offset = 0;
Def(1).SWIMS_num = [];
Idx(1).dtnum_beg = [];
Idx(1).dtnum_end = [];
Idx(1).filename = [];
PROG ='cruiseStruct_init_ADCP.m';

%cr_dir=pwd;
%DATA_dir = '/Users/johnmickett/Cruises_Research/PAPA/PAPA08/ADCP/Data/SN4021/data_mat';

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