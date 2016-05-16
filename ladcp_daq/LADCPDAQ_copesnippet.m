%For the next version of the LADCP software we'll keep the last x records
%of Vel.  So, whenever we have a new vel structure, we 

%Changes needed to GUI: 
%add std dev of btrange as well as mean
%add other beams to intensity plot in main gui
%add text widgets for zrange-min and zrange-max
%variables: handles.zrange_min, handles.zrange_max
%Make sure these update the plotting range when changed.

%Changed needed upon new run: 
%new variable, #ensembles in buffer - initialize to default value of 1000
%(about 1/2 hour). name: handles.plot_ensembles
%new variable, Vela, with this number of records, initialized to nan - wait
%till data come in

%open and clear plot window 1 - set position for PC screen

%needed files: 
%ezcf.m
%

%At each new bit of data:
%concatenate 
%    fn=fieldnames(str1);
%    for c=1:length(fn)
%        eval(['str3.' fn{c} '=[str1.' fn{c} ' str2.' fn{c} '];'])
%    end

%% Just look near the bottom approach
handles.plot_ensembles=1000;
fields={'depth_xducer';'yday';'dtnum';'heading';'pitch';'roll';'btrange';'ec1_bm';'u';'v';'btvel_bm';'z_adcp';'ens_no'};
%Do this the first time we have a valid Vel - this way we know the correct
%size.
clear Vela
    for c=1:length(fields)
        [m,n]=size(Vel.(fields{c}));
        Vela.(fields{c})=nan(m,handles.plot_ensembles);
    end
    
fields2={'vb';'vb_bt';'ub';'ub_bt';'ec1b'};
%Now make the depth vector
Vela.zg=(handles.zrange_min:4:handles.zrange_max)';
%Then make the fields.  This gets called when we start, or when we change
%the grid
    for c=1:length(fields2)
        Vela.(fields2{c})=nan(length(Vela.zg),handles.plot_ensembles);
    end

%%
%This little snippet should add the needed fields of Vel to Vela, removing
%the first record.
for c=1:length(fields)
    Vela.(fields{c})=[Vela.(fields{c})(:,2:end) Vel.(fields{c})];
end
%%
for c=1:length(fields)
    Vela.(fields{c})=Vel.(fields{c})(:,ig);
end

%%
%zrange=700;% meters above bottom to plot

%The approach for the LADCP-DAQ program should be to plot a specified dbar
%range; default to 4500-5200 m
%zrange = [4500 5200];

prdown_ff=1.88;

%ig=find((max(Vel.depth_xducer) - Vel.depth_xducer)*prdown_ff < zrange);
%ig=find((max(Vel.depth_xducer) - Vel.depth_xducer)*prdown_ff < zrange);
%For LADCP_DAQ, ig can be the whole thing
ig=1:length(Vela.yday);

%Then ig should find the last set time worth of Vel. Each ensemble will be
%added to the Vel structure.  To assure no memory problems, preallocate a
%structure with the needed fields and the alloted size.


%indicate that we need to regrid
handles.interp_start=1;


%interpolate each profile onto the grid accounting for the pressure of the
%ADCP

if handles.interp_start==length(ig)
end

for c=1:length(fields2)
    Vela.(fields2{c})=[Vela.(fields2{c})(:,2:end) nan(size(Vela.zg))];
end


for c=1:length(ig)
    Vela.ub(:,ig(c))=interp1(prdown_ff*Vela.depth_xducer(ig(c)) + Vela.z_adcp,Vela.u(:,ig(c)),Vela.zg);
    Vela.ub_bt(:,ig(c))=interp1(prdown_ff*Vela.depth_xducer(ig(c)) + Vela.z_adcp,Vela.u(:,ig(c))-Vela.btvel_bm(1,ig(c)),Vela.zg);
    Vela.vb(:,ig(c))=interp1(prdown_ff*Vela.depth_xducer(ig(c)) + Vela.z_adcp,Vela.v(:,ig(c)),Vela.zg);
    Vela.vb_bt(:,ig(c))=interp1(prdown_ff*Vela.depth_xducer(ig(c)) + Vela.z_adcp,Vela.v(:,ig(c))-Vela.btvel_bm(2,ig(c)),Vela.zg);
    Vela.ec1b(:,ig(c))=interp1(prdown_ff*Vela.depth_xducer(ig(c)) + Vela.z_adcp,Vela.ec1_bm(:,ig(c)),Vela.zg);
end

%indicate that we only need to do the next one
handles.interp_start=length(ig);
%% Plot - a dtnum version
figure(3)
clf
cl=[-.45 .45];
cl_i=[0 200];
ncolors=16;
xl=[min(Vela.dtnum(ig)) max(Vela.dtnum(ig))];
ax=MySubplot(.1,.1,0,.1,.1,0.02,1,8)
%
axes(ax(1))
plot(Vela.dtnum(ig),Vela.heading(ig));
xlim(xl)
datetick('x','keeplimits','keepticks')
xtloff
SubplotLetter('Heading',.01,.9)
ylim([0 360])
set(gca,'ytick',[0 90 180 270 360])
grid
axes(ax(2))
plot(Vela.dtnum(ig),Vela.pitch(ig), 'b-',Vela.dtnum(ig),Vela.roll(ig),'g-');
xlim(xl)
datetick('x','keeplimits','keepticks')
xtloff
SubplotLetter('pitch (blue), roll (green)',.01,.9)
grid
axes(ax(3))
%
ezcf(Vela.dtnum(ig),Vela.z_adcp,Vela.ec1_bm(:,ig),cl_i,ncolors);
%
hold on
plot(Vela.dtnum(ig),Vela.depth_xducer(ig)*prdown_ff);
hold off
%
%caxis([-.25 .25])
xlim(xl)
datetick('x','keeplimits','keepticks')
xtloff
%
ylabel('Range / m')
SubplotLetter('Beam 1 intensity',.01,.9)
axes(ax(4))

ezcf(Vela.dtnum(ig),Vela.zg,Vela.ec1b(:,ig),cl_i,ncolors)
hold on
plot(Vela.dtnum(ig),Vela.depth_xducer(ig)*prdown_ff,'k-')
plot(Vela.dtnum(ig),Vela.depth_xducer(ig)*prdown_ff+nanmean(Vela.btrange(:,ig))/100,'k-');
hold off
xlim(xl)
ylabel('P / dbar')
datetick('x','keeplimits','keepticks')
xtloff
SubplotLetter('Beam 1 intensity',.01,.9)
%
axes(ax(5))
ezcf(Vela.dtnum(ig),Vela.zg,Vela.ub(:,ig),cl,ncolors)
hold on
plot(Vela.dtnum(ig),Vela.depth_xducer(ig)*prdown_ff,'k-')
plot(Vela.dtnum(ig),Vela.depth_xducer(ig)*prdown_ff+nanmean(Vela.btrange(:,ig))/100,'k-');
xlim(xl)
datetick('x','keeplimits','keepticks')
xtloff
SubplotLetter('Measured u',.01,.9)
ylabel('P / dbar')
axes(ax(6))

ezcf(Vela.dtnum(ig),Vela.zg,Vela.vb(:,ig),cl,ncolors)
hold on
plot(Vela.dtnum(ig),Vela.depth_xducer(ig)*prdown_ff+nanmean(Vela.btrange(:,ig))/100,'k-');
xlim(xl)
datetick('x','keeplimits','keepticks')
SubplotLetter('Measured v',.01,.9)
ylabel('P / dbar')
xtloff

axes(ax(7))
ezcf(Vela.dtnum(ig),Vela.zg,Vela.ub_bt(:,ig),cl,ncolors)
hold on
plot(Vela.dtnum(ig),Vela.depth_xducer(ig)*prdown_ff,'k-')
plot(Vela.dtnum(ig),Vela.depth_xducer(ig)*prdown_ff+nanmean(Vela.btrange(:,ig))/100,'k-');
xlim(xl)
datetick('x','keeplimits','keepticks')
xtloff
SubplotLetter('u w/ bottom track',.01,.9)
ylabel('P / dbar')
axes(ax(8))

ezcf(Vela.dtnum(ig),Vela.zg,Vela.vb_bt(:,ig),cl,ncolors)
hold on
plot(Vela.dtnum(ig),Vela.depth_xducer(ig)*prdown_ff+nanmean(Vela.btrange(:,ig))/100,'k-');
xlim(xl)
datetick('x','keeplimits','keepticks')
SubplotLetter('v w/ bottom track',.01,.9)
ylabel('P / dbar')

linkaxes(ax,'x')

xlabel('Time UTC')
