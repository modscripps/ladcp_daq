function Vel = Run_WHbeams_mc09_AVG(yday_beg, yday_end, WHtyp, Params)
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
%           StartYear= starting year of deployment
%          NomDepth= nominal bottom depth at mooring location.


%here's where we can do some sort of for loop to average chunks of
%data...but would this be best to do in get_ADCP_Any? Must be sure to
%average onto a constant time grid that we overlap sections to avoid averaging problems near
%the ends of data chunks

%taking ten-minute boxcar averages of data and sub-sampling at 5 min
%intervals to avoid aliasing.
StartYear=2009; %add this to params.

Vel=[];
%breaking total into intervals of 1 day if longer than this.
TgridInt=5; %gridding interval in minutes
TMEgrid=[floor(yday_beg.*24)./24:TgridInt./1440:ceil(yday_end.*24)./24];
depthvec=[];

InTT=yday_end-yday_beg;
OVRLAP=(TgridInt.*4)./1440; %overlapping by 20 minutes so we can average and then join segments without edge problems
if InTT>1;  %if greater than a day break down interval into day-long segments for averaging
    IsLN=ceil(InTT); %number of sub-divisions
    for jj=1:IsLN;
        jj

        if jj==1; %if first sub-division use yday_beg
            yday_begS=yday_beg+(jj-1);
        else
            yday_begS=yday_beg+(jj-1-OVRLAP);
        end

        if jj==IsLN; %if last sub-division use yday_end
            yday_endS=yday_end;
        else
            yday_endS=yday_beg+jj+(OVRLAP);  %if not last sub-division increment by one day.
        end

        ADP = [];
        % WHtyp must be one of the strings in cell array Mnms:
        Mnms = {'SN3160m','SN3160a','SN11181m','SN7966m','SN10010m','SN11675m','SN8064m','SN7753m','SN11681m','SN4019m'};
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
        ADnms = {'SN3160','SN3160','SN11181','SN7966','SN10010','SN11675','SN8064','SN7753','SN11681','SN4019'};
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
                MooringID='WH1';
                Lat=36+(48.060./60);
                Lon=-[122+(48.648./60)];
                NomDepth=152.7;
            case 'SN3160a'
                ZGrid = [4:4:200]'; % Ascension Canyon (deploy#2)
                MagDec=14.1;
                MooringID='AC1'
                Lat=37+(01.468./60);
                Lon=-[122+(24.539./60)];
                NomDepth=199.25;
            case 'SN11181m'
                ZGrid = [4:8:324]'; % Long Ranger
                MooringID='LR4';
                Lat=[36+(47.676./60) 36+(48.032./60)];
                Lon=[-[121+(50.988./60)] -[121+(51.248./60)]];
                NomDepth=[307, 335.5];
                Notes='mooring was dragged 0.4 nm to NW on 03 Apr 09';
            case 'SN7966m'
                ZGrid = [2:2:50]';
                MooringID='MP1';
                Lat=36+(46.546./60);
                Lon=-[121+(55.4./60)];
                NomDepth=489;
            case 'SN10010m'
                ZGrid = [2:2:50]';
                MooringID='MP2';
                Lat=36+(47.232./60);
                Lon=-[121+(53.583./60)];
                NomDepth=369.9;
            case 'SN11675m'
                ZGrid = [2:2:50]';
                MooringID='MP3';
                Lat=36+(47.433./60);
                Lon=-[121+(52.038./60)];
                NomDepth=287.7;
            case 'SN8064m'
                %ZGrid = [32:4:220]'; % testing this unit (downward)
                ZGrid = [30:4:222]'; % testing this unit (downward)
                UpDown = 1;
                MooringID='MP1';
                Lat=36+(46.546./60);
                Lon=-[121+(55.4./60)];
                NomDepth=489;
                
               case 'SN7753m'
                ZGrid = [16:16:592]'; % testing this unit (downward)
                UpDown = 0;
                MooringID='LR1';
                Lat=36+(46.945./60);
                Lon=-[121+(55.115./60)];
                NomDepth=583.8;
                
                  case 'SN11681m'
                ZGrid = [16:16:578]'; % testing this unit (downward)
                UpDown = 0;
                MooringID='LR2';
                Lat=36+(46.416./60);
                Lon=-[121+(54.007./60)];
                NomDepth=578.9;
                
                  case 'SN4019m'
                ZGrid = [8:8:408]'; % testing this unit (downward)
                UpDown = 0;
                MooringID='LR3';
                Lat=36+(47.370./60);
                Lon=-[121+(53.994./60)];
                NomDepth=419.8;             
            otherwise
                disp(['Run_WHbeams_mc09:  WHtyp=' WHtyp ' is unkwown, exit'])
                return
        end
        if exist('Params', 'var') && isstruct(Params)
            vvs = {'MatFld','MatIndx','WC_val','Psrc','Ztyp','ZGrid','Sal0','MagDec','NomDepth','StartYear'};
            for i=1:length(vvs)
                if isfield(Params, vvs{i}) && ~isempty(Params.(vvs{i}))
                    eval([vvs{i} '=Params.(vvs{i});']);
                end
            end
        end
        %
        % Now, get beam velocities for this ADCP
        ADRaw = get_ADCP_Any(yday_begS, yday_endS, MatIndx, MatFld, 1);


        if isempty(ADRaw)
            disp(['Run_WHbeams_mc09: Cannot get ' WHtyp ' data, exit'])
            continue
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
            %ADP=rmfield(ADP,'z_adcp');

        else % stg1: no soundspeed correction, ADCP-relative depth bins
            WHbeam_Process_stg1
        end

        %now we're averaging the ADP data

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
        PInt=round(nanmedian(diff(ADP.yday.*86400).*10))./10;

        %culling data that exceeds werr threshold
        Isbb=find(abs(ADP.werr>1));
        ADP.u(Isbb)=NaN;
        ADP.v(Isbb)=NaN;
        ADP.w(Isbb)=NaN;
        
        
        
        MagVel=sqrt(ADP.u.^2+ADP.v.^2);
        MagVelR=reshape(MagVel,1,size(MagVel,1).*size(MagVel,2));
        Bstd=nanstd(MagVelR);
      
        
%          Ibdd=find(MagVel>=(nanmean(MagVelR)+3.*Bstd));
%          
%          ADP.u(Ibdd)=NaN;
%          ADP.v(Ibdd)=NaN;
%          ADP.w(Ibdd)=NaN;
%          clear Ibdd
%         
%       Ibdd=find(MagVel>=0.7);
%          
%          ADP.u(Ibdd)=NaN;
%          ADP.v(Ibdd)=NaN;
%          ADP.w(Ibdd)=NaN;   
        
        

        ecMM=ones(length(ADP.depth),length(ADP.yday),4);
        %averaging echo intensities
        ecMM(:,:,1)=ADP.ec1;
        ecMM(:,:,2)=ADP.ec2;
        ecMM(:,:,3)=ADP.ec3;
        ecMM(:,:,4)=ADP.ec4;
        ecM=nanmean(ecMM,3);


        %getting rid of data w/in 10% of sampling depth or near bottom
        %taking 10% of depth xducer, removing this from top

        if nanmean(ADP.z_adcp<0); %upward-looking
            Xblank=0.1.*nanmax(ADP.depth_xducer);
            Ido=find(ADP.depth<=Xblank);
        else %downward looking
            %cutting off velocities below bottom and those within 10% of bottom
            DiffD=NomDepth-nanmedian(ADP.depth_xducer);
            GoodD=0.9.*DiffD;
            Ido=find(ADP.z_adcp>(nanmedian(ADP.depth_xducer)+GoodD));
        end

        ADP.u(Ido,:)=NaN;
        ADP.v(Ido,:)=NaN;
        ADP.w(Ido,:)=NaN;
        
  
        if ~isempty(ADRaw) & isempty(depthvec);
        depthvec=ADP.depth;
        end
        




        %if it is the first day, creating matrices
        if ~exist('uMat')==1 & ~isempty(ADP);

            uMat=NaN.*ones(length(ADP.depth),length(TMEgrid));
            vMat=uMat;
            wMat=uMat;
            ecMat=uMat;
            xDucerDep=ones(1,length(TMEgrid));
            tempM=xDucerDep;

        end






        %creating 10 minute boxcar regardless of sampling interval

        SampInt=nanmean(diff(ADP.dtnum)).*86400; %sampling interval in seconds
        NNB=ceil((10.*60)./SampInt); %rounding up...this is the number of samples in the 10-minute average.
        BB=boxcar(NNB)./nansum(boxcar(NNB));


        %averaging variables here



        uAve=conv2(1,BB',ADP.u,'same');
        vAve=conv2(1,BB',ADP.v,'same');
        wAve=conv2(1,BB',ADP.w,'same');
        ecMAve=conv2(1,BB',ecM,'same');
        xducer_depth=conv2(1,BB',ADP.depth_xducer,'same');
        tempAve=conv2(1,BB',ADP.temp,'same');

        %NaNing beginning and ends where averaging causes problems. Being overly
        %conservative here.
        uAve(:,1:NNB)=NaN;
        uAve(:,end-NNB:end)=NaN;
        vAve(:,1:NNB)=NaN;
        vAve(:,end-NNB:end)=NaN;
        wAve(:,1:NNB)=NaN;
        wAve(:,end-NNB:end)=NaN;
        ecMAve(:,1:NNB)=NaN;
        ecMAve(:,end-NNB:end)=NaN;

        xducer_depth(end-NNB:end)=NaN;
        xducer_depth(1:NNB)=NaN;
        tempAve(end-NNB:end)=NaN;
        tempAve(1:NNB)=NaN;


        Isff=find(TMEgrid>=ADP.yday(1) & TMEgrid<=ADP.yday(end)); %to save time only interpolating onto a
        %portion of the total time.

        %interpolating onto standard grid.

        for ii=1:length(ADP.depth);
            %ii=10;
            %u
            [C,Iuu]=unique(ADP.yday);

            uAveI=interp1(ADP.yday(Iuu),uAve(ii,Iuu),TMEgrid(Isff),'linear'); %interpolating onto 5-min grid after smoothing
            iso=find(~isnan(uAveI));
            uMat(ii,Isff(iso))=uAveI(iso);
            %v
            vAveI=interp1(ADP.yday(Iuu),vAve(ii,Iuu),TMEgrid(Isff),'linear'); %interpolating onto 5-min grid after smoothing
            iso=find(~isnan(vAveI));
            vMat(ii,Isff(iso))=vAveI(iso);
            %w
            wAveI=interp1(ADP.yday(Iuu),wAve(ii,Iuu),TMEgrid(Isff),'linear'); %interpolating onto 5-min grid after smoothing
            iso=find(~isnan(wAveI));
            wMat(ii,Isff(iso))=wAveI(iso);
            %ec
            ecAveI=interp1(ADP.yday(Iuu),ecMAve(ii,Iuu),TMEgrid(Isff),'linear'); %interpolating onto 5-min grid after smoothing
            iso=find(~isnan(ecAveI));
            ecMat(ii,Isff(iso))=ecAveI(iso);

        end

        xducerDepthI=interp1(ADP.yday(Iuu),xducer_depth(Iuu),TMEgrid(Isff),'linear'); %interpolating onto 5-min grid after smoothing
        xDucerDep(Isff(iso))=xducerDepthI(iso);


        tempAveI=interp1(ADP.yday(Iuu),tempAve(Iuu),TMEgrid(Isff),'linear'); %interpolating onto 5-min grid after smoothing
        tempM(Isff(iso))=tempAveI(iso);



        %putting pieces back together, remembering we overlapped by 20 minutes, so just
        %adding unique time values.

    end
end

%finding all nans at beginning or end of timeseries, removing from
%timeseries

%this gets all nans, just removing beginning and end sections
Igood=find(~isnan(nanmean(uMat)));

if ~isempty(Igood);
IgdS=Igood(1);
IgdE=Igood(end);
else
    disp('No good data within range specified');
    return
end


%if LR4 then create a vector of lat and lon AND of nominal water depth.
%with depth of Xducer 11 m above the bottom.





Vel.u=uMat(:,IgdS:IgdE);
Vel.v=vMat(:,IgdS:IgdE);
Vel.w=wMat(:,IgdS:IgdE);
Vel.ecAve=ecMat(:,IgdS:IgdE);
Vel.yday=TMEgrid(IgdS:IgdE);
Vel.datenum=yday2datenum(TMEgrid(IgdS:IgdE),StartYear);
Vel.z=depthvec;
Vel.xducer_depth=xDucerDep(IgdS:IgdE);
Vel.H=NomDepth;
Vel.lat=Lat;
Vel.lon=Lon;
Vel.temp=tempM(IgdS:IgdE);
Vel.info.SNadcp=WHtyp(1:end-1);
Vel.info.station=MooringID;
Vel.info.cruise='MC09';
Vel.info.startyear=StartYear;
Vel.info.mfiles=['in order: Conv_ADCP_mooring (coverts to mat files)'...
    'Run_WHbeams_mc09_AVG.m (conv. to earth coord +avgs. & calls get_ADCP_Any.m, WHbeam_Process.m, etc.)'];
Vel.info.processing=['10-minute boxcar avgs. from ' num2str(PInt) 's pings, subsampled at ' num2str(TgridInt) ' minutes.'];
if exist('Notes')==1;
    Vel.info.notes=Notes;
end

if strcmp(WHtyp,'SN11181m')==1; %do long ranger stuff.
II=find(Vel.xducer_depth>310);
NewLat=0.*Vel.yday+Vel.lat(1);
NewLon=0.*Vel.yday+Vel.lon(1);
NewH=0.*Vel.yday+Vel.H(1);
if ~isempty(II);
NewLat(II)=Vel.lat(2).*ones(1,length(II));
NewLon(II)=Vel.lon(2).*ones(1,length(II));
NewH(II)=Vel.H(2).*ones(1,length(II));
end
Vel.H=NewH;
Vel.lon=NewLon;
Vel.lat=NewLat;
end

