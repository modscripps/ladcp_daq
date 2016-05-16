%this script merges two ADCP records from the same mooring
%data should be fully processed first for each ADCP.


%loading data:


%from PAPA08 mooring.
load('/Volumes/Puao/data_archive/PAPA/PAPA08/ADCP/Data/SN0578/data_mat/SN0578_PAPA08_All.mat');
Vel1=Vel;
load('/Volumes/Puao/data_archive/PAPA/PAPA08/ADCP/Data/SN4021/data_mat/SN4021_PAPA08_All.mat');
Vel2=Vel;
clear Vel


% Vel2 = 
%                u: [51x17705 double]
%                v: [51x17705 double]
%                w: [51x17705 double]
%            ecAve: [51x17705 double]
%             yday: [1x17705 double]
%            dtnum: [1x17705 double]
%                z: [51x1 double]
%     xducer_depth: [1x17705 double]
%             temp: [1x17705 double]
%         botdepth: 4235
%              lat: 50
%              lon: -145
%             info: [1x1 struct]
%setting time and depth grid.

%gridding to 
Iee=(nanmedian(diff(Vel1.dtnum))); %hourly on min 53
Igg=(nanmedian(diff(Vel2.dtnum))); %half-hourly on 1/2 hours
%gridding to finer grid so as not to lose data

Dtvec=Vel2.dtnum;

%preserving vertical resolution, so patching both z-grids.


Ise=find(Vel2.z>Vel1.z(end));

U2new=Vel2.u(Ise,:);
Z2new=Vel2.z(Ise);
V2new=Vel2.v(Ise,:);
W2new=Vel2.w(Ise,:);
Ec2new=Vel2.ecAve(Ise,:);

Znew=[Vel1.z;Z2new];

V1I=interp2(Vel1.dtnum,Vel1.z,Vel1.v,Dtvec,Vel1.z);
U1I=interp2(Vel1.dtnum,Vel1.z,Vel1.u,Dtvec,Vel1.z);
W1I=interp2(Vel1.dtnum,Vel1.z,Vel1.w,Dtvec,Vel1.z);
Ec1I=interp2(Vel1.dtnum,Vel1.z,Vel1.ecAve,Dtvec,Vel1.z);


Vmrg=[V1I;V2new];
Umrg=[U1I;U2new];
Wmrg=[W1I;W2new];
Ecmrg=[Ec1I;Ec2new];


figure
pcolorjbm(Dtvec,Znew,Vmrg);caxis([-0.5 0.5]);colormap(redblue);


figure
pcolorjbm(Dtvec,Znew,Ecmrg);caxis([0 255]);colormap(jet);
%OK, now putting onto a uniform depth grid (4 m)

ZgridN=[Znew(1):4:Znew(end)]'; %this preserves vertical resolution in top 150 m.

Utot=NaN.*ones(length(ZgridN),length(Dtvec));
Vtot=Utot;
Wtot=Utot;
Ectot=Utot;


Umrga=Umrg;
Vmrga=Vmrg;
Wmrga=Wmrg;
Ecmrga=Ecmrg;

Umrg(38:41,:)=NaN;
Vmrg(38:41,:)=NaN;
Wmrg(38:41,:)=NaN;
Ecmrg(38:41,:)=NaN;


%linearly interpolating over NaN gaps
for ii=1:length(Dtvec);
    Ibc=find(~isnan(Umrg(:,ii)));
    if length(Ibc)>10;
    uta=interp1(Znew(Ibc),Umrg(Ibc,ii),ZgridN,'linear');
    Utot(:,ii)=uta;
    end
    
    Ibd=find(~isnan(Vmrg(:,ii)));
    if length(Ibd)>10;
    vta=interp1(Znew(Ibd),Vmrg(Ibd,ii),ZgridN,'linear');
    Vtot(:,ii)=vta;
    end
    
    Ibe=find(~isnan(Wmrg(:,ii)));
    if length(Ibe)>10;
    wta=interp1(Znew(Ibe),Wmrg(Ibe,ii),ZgridN,'linear');
    Wtot(:,ii)=wta;
    end
    
      Ibf=find(~isnan(Ecmrg(:,ii)));
    if length(Ibf)>10;
    ecta=interp1(Znew(Ibf),Ecmrg(Ibf,ii),ZgridN,'linear');
    Ectot(:,ii)=ecta;
    end
    
    
    
    
end


temp1=interp1(Vel1.dtnum,Vel1.temp,Dtvec);
xdepth1=interp1(Vel1.dtnum,Vel1.xducer_depth,Dtvec);


% figure
% pcolorjbm(Dtvec,ZgridN,Vtot);caxis([-0.5 0.5]);colormap(redblue);


% Vel2 = 
%                u: [51x17705 double]
%                v: [51x17705 double]
%                w: [51x17705 double]
%            ecAve: [51x17705 double]
%             yday: [1x17705 double]
%            dtnum: [1x17705 double]
%                z: [51x1 double]
%     xducer_depth: [1x17705 double]
%             temp: [1x17705 double]
%         botdepth: 4235
%              lat: 50
%              lon: -145
%             info: [1x1 struct]
%setting time and depth grid.




Vel.u=Utot;
Vel.v=Vtot;
Vel.w=Wtot;
Vel.ecAve=Ectot;

Vel.yday=Vel2.yday;
Vel.dtnum=Dtvec;
Vel.z=ZgridN;
Vel.temps=[temp1;Vel2.temp];
Vel.xducer_depths=[xdepth1;Vel2.xducer_depth];
Vel.botdepth=Vel1.botdepth;
Vel.lat=Vel1.lat;
Vel.lon=Vel1.lon;
Vel.info=Vel1.info;
Vel.info.sysconfig.freq=[{Vel1.info.sysconfig.freq};{Vel2.info.sysconfig.freq}];;
Vel.info.sysconfig.beamangle=[];
Vel.info.snADCP=[{Vel1.info.snADCP};{Vel2.info.snADCP}];
Vel.info.NomInstDepth=[{Vel1.info.NomInstDepth};{Vel2.info.NomInstDepth}];
Vel.info.mfiles='merger_moored_ADCPs.m';

VeluI=Vel;
VeluI.u=Umrga;
VeluI.v=Vmrga;
VeluI.w=Wmrga;
VeluI.z=Znew;


