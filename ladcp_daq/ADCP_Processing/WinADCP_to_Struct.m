function Vel = WinADCP_to_Struct(in_file,Lat,depth_xducer);
%INPUT in_file is the mat file output from WinADCP, depth_xducer is the depth of the transducer in m (a single value)
%or it is a 2x1 vector of [dpth dtnum].,
%Lat is latitude in deg. for computing pressure from depth using sw_pres.m
%this loads the mat file exported from WinADCP and saves it to a standard
%ADCP structure for raw-ish data-, for earth coords, declination NOT
%corrected for, etc.)  JBM 7/1/09...

WV=load(in_file);


% %WV = 
%        RDIFileName: [1x55 char]
%          RDISystem: 'Narrowband 307.2 kHz'
%     RDIPingsPerEns: 300
%      RDISecPerPing: 3
%         RDIEnsDate: '05/30'
%         RDIEnsTime: '23:53:00.00'
%     RDIEnsInterval: 1200
%         RDIBin1Mid: 6
%         RDIBinSize: 4
%            SerBins: [1x44 double]
%       SerEnsembles: [9166x1 double]
%            SerYear: [9166x1 double]
%             SerMon: [9166x1 double]
%             SerDay: [9166x1 double]
%            SerHour: [9166x1 double]
%             SerMin: [9166x1 double]
%             SerSec: [9166x1 double]
%            SerHund: [9166x1 double]
%        AnP100thDeg: [9166x1 double]
%        AnR100thDeg: [9166x1 double]
%        AnH100thDeg: [9166x1 double]
%        AnT100thDeg: [9166x1 double]
%          AnDepthmm: [9166x1 double]
%          AnOrienUP: [9166x1 double]
%              AnBIT: [9166x1 double]
%             AnBatt: [9166x1 double]
%          SerEA1cnt: [9166x44 double]
%          SerEA2cnt: [9166x44 double]
%          SerEA3cnt: [9166x44 double]
%          SerEA4cnt: [9166x44 double]
%          SerEAAcnt: [9166x44 double]
%           SerC1cnt: [9166x44 double]
%           SerC2cnt: [9166x44 double]
%           SerC3cnt: [9166x44 double]
%           SerC4cnt: [9166x44 double]
%           SerCAcnt: [9166x44 double]
%       SerEmmpersec: [9166x44 double]
%       SerNmmpersec: [9166x44 double]
%       SerVmmpersec: [9166x44 double]
%      SerErmmpersec: [9166x44 double]
%     SerMagmmpersec: [9166x44 double]
%      SerDir10thDeg: [9166x44 double]
%             SerPG1: [9166x44 double]
%             SerPG2: [9166x44 double]
%             SerPG3: [9166x44 double]
%             SerPG4: [9166x44 double]
% 
% 
% % Vel = 
% %           ens_no: [1x18093 double]
% %           z_adcp: [1x55 double]
% %           p_adcp: [1x55 double]
% %         pulselen: 1673
% %             yday: [1x18093 double]


% %          heading: [1x18093 double]
% %            pitch: [1x18093 double]
% %             roll: [1x18093 double]
% %          hdg_std: [1x18093 double]
% %        pitch_std: [1x18093 double]
% %         roll_std: [1x18093 double]
% %     depth_xducer: [1x18093 double]
% %         soundvel: [1x18093 double]
% %             temp: [1x18093 double]
% %               v1: [55x18093 double]
% %               v2: [55x18093 double]
% %               v3: [55x18093 double]
% %               v4: [55x18093 double]
% %           ec1_bm: [55x18093 double]
% %           ec2_bm: [55x18093 double]
% %           ec3_bm: [55x18093 double]
% %           ec4_bm: [55x18093 double]
% %          cor1_bm: [55x18093 double]
% %          cor2_bm: [55x18093 double]
% %          cor3_bm: [55x18093 double]
% %          cor4_bm: [55x18093 double]
% %              pg1: [55x18093 double]
% %              pg2: [55x18093 double]
% %              pg3: [55x18093 double]
% %              pg4: [55x18093 double]

%correcting strange problem with SerYear for NOAA 300 kHz narrowband (all
%values are 80!).
ID=99; %99 if PAPA08_0578 mooring, else something else

if ID==99;
    
Ido=find(abs(diff(WV.SerMon))==nanmax(abs(diff(WV.SerMon))));

SerYear=NaN.*WV.SerMon;

SerYear(1:Ido)=08;
SerYear(Ido+1:end)=09;
end


Vel.ens_no=WV.SerEnsembles;

SerYear=SerYear+2000;

WV.SerYear=SerYear;

dnum=datenum(WV.SerYear',WV.SerMon',(WV.SerDay)',WV.SerHour',WV.SerMin',(WV.SerSec+WV.SerHund./100)');





if WV.AnOrienUP(1)==1; %uplooking
Vel.z_adcp=-(WV.RDIBin1Mid+(WV.RDIBinSize.*[0:1:(length(WV.SerBins)-1)]));
else
Vel.z_adcp=WV.RDIBin1Mid+(WV.RDIBinSize.*[0:1:(length(WV.SerBins)-1)]);
end



Vel.p_adcp=[];
Vel.pulselen=[];

ydayS=yearday(WV.SerDay',WV.SerMon',(WV.SerYear)',WV.SerHour',WV.SerMin',(WV.SerSec+WV.SerHund./100)');

%Vel.yday=ydayS+(dnum-dnum(1));


Vel.dtnum=dnum;

Vel.yday=datenum2yday(dnum);




Vel.heading=WV.AnH100thDeg'./100;
Vel.pitch=WV.AnP100thDeg'./100;
Vel.roll=WV.AnR100thDeg'./100;
Vel.heading_std=[];
Vel.pitch_std=[];
Vel.roll_std=[];

if nargin>2 & length(depth_xducer)>1;  %if there is a vector of depth inputs

    dxducer=interp1(depth_xducer(:,2),depth_xducer(:,1),Vel.dtnum);
    
Vel.depth_xducer=dxducer;
elseif nargin>2 & length(depth_xducer)==1;
    dxducer=NaN.*Vel.dtnum+depth_xducer;
    Vel.depth_xducer=dxducer;
else
    Vel.depth_xducer=WV.AnDepthmm./1000;
end

Sal=35;

PP= sw_pres(nanmean(Vel.depth_xducer),Lat);

Vel.soundvel = sw_svel(Sal./1000,nanmedian(WV.AnT100thDeg./100),nanmedian(PP)./100);

Vel.temp=WV.AnT100thDeg'./100;



Vel.v1=WV.SerEmmpersec'./1000;
Vel.v2=WV.SerNmmpersec'./1000;

Vel.v3=WV.SerVmmpersec'./1000;
Vel.v4=WV.SerErmmpersec'./1000; %is this the same?

Vel.ec1_bm=WV.SerEA1cnt';
Vel.ec2_bm=WV.SerEA2cnt';
Vel.ec3_bm=WV.SerEA3cnt';
Vel.ec4_bm=WV.SerEA4cnt';

Vel.cor1_bm=WV.SerC1cnt';
Vel.cor2_bm=WV.SerC2cnt';
Vel.cor3_bm=WV.SerC3cnt';
Vel.cor4_bm=WV.SerC4cnt';

Vel.pg1=WV.SerPG1';
Vel.pg2=WV.SerPG2';
Vel.pg3=WV.SerPG3';
Vel.pg4=WV.SerPG4';


Vel.sysconfig.beamdir='upfacing';
Vel.sysconfig.beampattern='convex';
Vel.sysconfig.freq='300kHz';
Vel.sysconfig.beamangle='20?';
Vel.params=[];
Vel.mfiles=[];