%LADCP_DAQ.m
%Version 1.1, MHA
%Written to work with MATLAB r2008a on a windows machine to use the serial
%toolbox to monitor a serial stream and write it to disk.
%We use GUIDE to make a GUI window to control and view data aquisition.
%
%We use the serial toolbox to monitor the com1 port for incoming hex data
%from the LADCP as output from the downlooking ADCP via the SBE 9 uplink.
%Whenever we receive data, we search it for beginnings of RDI ensemble
%records and write them one by one to a binary file.  The base filenum is
%specified by the user, the extension .000 etc is incremented when each
%file reaches a certain size.  After each valid record is written, it is
%decoded using Winkelware and basic values and plots are made.

%v1.1 has some additional plotting.
%MHA 8/8/2012
%
%MHA 1/28/2014: added bottom warnings.

function varargout = LADCP_DAQ(varargin)
% LADCP_DAQ M-file for LADCP_DAQ.fig
%      LADCP_DAQ, by itself, creates a new LADCP_DAQ or raises the existing
%      singleton*.
%
%      H = LADCP_DAQ returns the handle to a new LADCP_DAQ or the handle to
%      the existing singleton*.
%
%      LADCP_DAQ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LADCP_DAQ.M with the given input arguments.
%
%      LADCP_DAQ('Property','Value',...) creates a new LADCP_DAQ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LADCP_DAQ_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LADCP_DAQ_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LADCP_DAQ

% Last Modified by GUIDE v2.5 10-Aug-2012 00:51:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LADCP_DAQ_OpeningFcn, ...
    'gui_OutputFcn',  @LADCP_DAQ_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before LADCP_DAQ is made visible.
function LADCP_DAQ_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LADCP_DAQ (see VARARGIN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MHA: begin MHA initialization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set up the serial port
handles.serial_port = 'com1'; %Can be set to a USB port as well.
handles.serial_s = serial(handles.serial_port);
handles.serial_s.BaudRate = 19200;
handles.serial_s.Terminator = '';
handles.serial_s.InputBufferSize=512; %default = 512.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial parameters

%Default output directory and file name.
%handles.outputdir='C:\Documents and Settings\All Users\Desktop\LADCP_DAQ\data_km1228\';
handles.outputdir='C:\Documents and Settings\All Users\Desktop\LADCP_DAQ\data_tgt305\';
%handles.outfilename='RR1209_LADCP_Downlooker_cast1';
handles.outfilename='TGT305_LADCP_Downlooker_cast1';

%initialize counters etc
handles.bufferbytes=0; %bytes in serial buffer
handles.filebytes=0; %bytes written to file
handles.lastline='';
handles.status='Idle.';

handles.filenum=0; %increment this after a certain number of records
handles.recsperfile=100; %number of records to store per file
handles.MAXFILEBYTES=5e6; %size of each file in bytes. 5 MB
handles.plotskip=3; %plot every (this)'th record; can set to greater than 1 if processing can't keep up

%plotting parameters
handles.zmin=3000;
handles.zmax=3500;
handles.start_ensemble=1;
handles.vmin=-0.2;
handles.vmax=0.2;
handles.enable_plots=0;

%enable bottom warnings
handles.enable_bottom_warnings=1;

%set GUI widgets to initial values
set(handles.outfilename_edit,'String',handles.outfilename);
set(handles.outputdir_edit,'String',handles.outputdir);
set(handles.lastline_text,'String',handles.lastline);
set(handles.status_text,'String',handles.status);

set(handles.bufferbytes_text,'String',num2str(handles.bufferbytes));
set(handles.filebytes_text,'String',num2str(handles.filebytes));

set(handles.edit_zmin,'String',num2str(handles.zmin));
set(handles.edit_zmax,'String',num2str(handles.zmax));

set(handles.edit_vmin,'String',num2str(handles.vmin));
set(handles.edit_vmax,'String',num2str(handles.vmax));

%set button enable status - start initially enabled.
set(handles.start_button,'Enable','On');
set(handles.stop_button,'Enable','Off');

%Add the path to the ADCP processing directory
addpath('ADCP_Processing')

% MHA: end MHA initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Choose default command line output for LADCP_DAQ
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LADCP_DAQ wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Executes on button press in start_button.
function start_button_Callback(hObject, eventdata, handles)
% hObject    handle to start_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Here we begin data aquisition by opening an output file and the serial port.

%First reset counters.
handles.bufferbytes=0; %bytes in serial buffer
handles.filebytes=0; %bytes written to file
handles.filenum=0;

%set up a bigger buffer than 512 bytes as RDI records are longer.
handles.buffa=[];
handles.count=0;

%Open output file.
handles.outfilename=get(handles.outfilename_edit, 'String');
%Check for existing files.
if 1 %SET TO 0 for testing
    if exist(fullfile(handles.outputdir, [handles.outfilename '.' sprintf('%03d',handles.filenum)]))==2
        set(handles.status_text,'String','File exists.  Please move or remove the file manually.')
        return;
    end
end

%Open file for writing; return if failure.
%Open it as a binary file with littleendian format as required by RDI.
handles.out = fopen(fullfile(handles.outputdir, [handles.outfilename '.' sprintf('%03d',handles.filenum)]),...
    'wb','l');

if handles.out < 0
    set(handles.status_text,'String','Error opening output file.')
    return;
end

%If we are here, we have a valid open file and are ready to open the serial
%stream and collect data.

%Set GUI's
set(handles.bufferbytes_text,'String',num2str(handles.bufferbytes));
set(handles.filebytes_text,'String',num2str(handles.filebytes));
set(handles.filenum_text,'String',num2str(handles.filenum));

%initialize GUI's in case they have info from last run in them
set(handles.lastline_text,'String','');
set(handles.ensemble_text,'String','');
set(handles.pressure_text,'String','');
set(handles.pitch_text,'String','');
set(handles.roll_text,'String','');
set(handles.heading_text,'String','');
set(handles.range_text,'String','');

%1/28/2014: initialize bottom detector
%intitialize
handles.bthits=0;
handles.btrange=nan;
handles.decreasing=0;

handles.max_warn_freq_sec= 60*4; %max frequency in secs of warnings.
handles.hits_to_warn=6; %required successive bottom detections to warn

handles.last=now - handles.max_warn_freq_sec/3600/24; %set the time of the last bottom warning in the past




%disable outfilename edit box so user cannot edit while we are aquiring
%data
set(handles.outfilename_edit,'Enable','Off');
set(handles.outputdir_edit,'Enable','Off');
set(handles.status_text,'String','Aquiring Data.')

%set button enable status
set(handles.start_button,'Enable','Off');
set(handles.stop_button,'Enable','On');

%Set stream properties.
handles.serial_s.BytesAvailableFcnCount = 512;
handles.serial_s.BytesAvailableFcnMode = 'byte';
handles.serial_s.BytesAvailableFcn = {@ParseData, hObject};

%Plotting parameters
handles.start_ensemble=1; %the first ensemble that needs to be gridded
handles.no_ensembles=300; %number of ensembles to store
%remove Vela if it exists
if isfield(handles,'Vela')
    handles=rmfield(handles,'Vela') %
end

%Open serial stream.
fopen(handles.serial_s);

% Update handles structure
guidata(hObject, handles);

function [] = ParseData(s, eventdata, fig)
%This is the Callback function when 512 bytes are received in the serial
%buffer.  We read the data, look for the RDI record identifier, write the
%record, and plot it.

%Get handles from the fig object
handles = guidata(fig);

%read synchronously from the serial port
buff = fscanf(s);

%add buff to the larger buffer
handles.buffa=[handles.buffa buff];

%update the number of bytes and its GUI
handles.bufferbytes=handles.bufferbytes+length(buff);
set(handles.bufferbytes_text,'String',num2str(handles.bufferbytes));

%Find the beginning of an RDI record.
%i1=findstr(handles.buffa,'7F7F0'); %8/12/12 mha: This string is for 8-m bins
%i1=findstr(handles.buffa,'7F7F93');  %16-m bins
%i1=findstr(handles.buffa,'7F7F47'); %1/11/13 mha: Testing after downlooker strangeness
i1=findstr(handles.buffa,'7F7F7F'); %1/14/13 mha: 20 bins, 16-m



if length(i1) >= 2 % do nothing unless we have a complete record, which we know from the presence of two 7F7F93 strings
    %put the data in the outbox
    outstr=handles.buffa(i1(1):i1(2)-1);
    
    %Make sure it is even length
    if rem(length(outstr),2)==1
        outstr=outstr(1:end-1);
    end
    %We need to turn the hex strings into binary bytes.  Each byte is
    %represented by two hex characters.  We'll do this by reshaping the matrix
    %into two rows and calling hex2dec on each column.
    
    %Then reshape so that each two hex characters are in a column (ie the first
    %and second from the original string are in the first column, the third and
    %fourth characters in the second column, and so on).
    tmp=reshape(outstr,2,length(outstr)/2)';
    
    %2/2/2014: check for mod errors which appear as non-hex chars
    mod_err_check=1;
    if mod_err_check
        moderr=any(any(~((tmp>='0' & tmp<='9') | (tmp>='A'&tmp<='F'))));
    else
        moderr=0;
    end
    if moderr
        disp 'mod error!'
            %set GUI with the record written
            set(handles.lastline_text,'String','');
            set(handles.lastline_text,'String',outstr);
            %remove the part written from the buffer
            handles.buffa=handles.buffa(i1(2):end);
            %update the buffer count
            handles.bufferbytes=handles.bufferbytes - length((1:i1(2)-1));

        moderr=0; %reset
    else %~moderr
        %Then turn this into an unsigned integer by calling hex2dec on each column
        %if any(any(~((tmp>='0' & tmp<='9') | (tmp>='A'&h<='F'))))
        %   error(message('MATLAB:hex2dec:IllegalHexadecimal'));
        %end
        out=hex2dec(tmp)';
        
        %Write to the file, return number of bytes written
        count=fwrite(handles.out,out,'uint8');
        
        %increment counter. Used to check if it is time to plot if we don't do
        %it every time
        handles.count=handles.count+1;
        handles.filebytes=handles.filebytes+count;
        set(handles.filebytes_text,'String',num2str(handles.filebytes));
        %    set(handles.filebytes_text,'String',num2str(count));
        
        %Output diagnostics to command window
        if 0
            disp(['#' num2str(handles.count) ':'])
            disp(['    length of buffa is ' num2str(length(handles.buffa))])
            %    disp(['    length of bufftmp is ' num2str(length(handles.bufftmp))])
            disp(['    length of buff is ' num2str(length(buff))])
            disp(['valid starts found in buffa at ' num2str(i1(1)) '&' num2str(i1(2)) '(diff=' num2str(diff(i1)) ')'])
            disp(['    length of output data is ' num2str(length(outstr))])
            disp(['    bytes of converted output data is ' num2str(length(out))])
            disp(['    bytes actually written is ' num2str(count)])
            disp(['    First bit to be written is ' outstr(1:10)])
            disp([num2str(length((1:i1(2)-1))) ' values removed from buffa'])
        end
        
        if count == length(out) %if a successful write
            %set GUI with the record written
            set(handles.lastline_text,'String','');
            set(handles.lastline_text,'String',outstr);
            %remove the part written from the buffer
            handles.buffa=handles.buffa(i1(2):end);
            %update the buffer count
            handles.bufferbytes=handles.bufferbytes - length((1:i1(2)-1));
            %This should be the same as length(handles.buffa); can check in
            %debugger
        else
            set(handles.lastline_text,'String','WARNING: less than expected bytes written.')
        end
        
        %Check to see if it's time to open a new output file
        if handles.filebytes > handles.MAXFILEBYTES
            %if so, close current file and open the new one
            fclose(handles.out);
            %increment the file number and update gui
            handles.filenum=handles.filenum+1;
            set(handles.filenum_text,'String',num2str(handles.filenum));
            
            %reset the bytes counter and update gui
            handles.filebytes=0;
            set(handles.filebytes_text,'String',num2str(handles.filebytes));
            
            %Open the new file.
            handles.out = fopen(fullfile(handles.outputdir, [handles.outfilename '.' sprintf('%03d',handles.filenum)]),...
                'wb','l');
            %if an error, try to exit gracefully.
            if handles.out < 0
                set(handles.status_text,'String','Error opening output file.')
                stop_button_Callback(fig, eventdata, handles)
            end
        end
        
        % Now that the output file is written (our first priority),
        % we can decode the last ensemble and make some plots etc. We
        % do this by writing a temporary file with a single ensemble, and then
        % using Winkel-ware to decode it into a Vel structure.
        if rem(handles.count,handles.plotskip)==0 %do every nth record (try 1; make sure we keep up with data...)
            %open a temporary file
            tmpfil = fopen('tmpfil.000','wb','l');
            %Write to the file, return number of bytes written
            count=fwrite(tmpfil,out,'uint8');
            %close file
            fclose(tmpfil);
            %Translate the file
            Vel=Get_ADCP_fullSC_BT('tmpfil.000',-1,1);
            %John's "new" version for LADCP
            %        Vel=Get_ADCP_fullSC_LADCP('tmpfil.000',-1,1);
            if ~isempty(Vel)
                %1/15/2014: first fix pressure wrap
                ibb=find(Vel.depth_xducer<0);
                if ~isempty(ibb)
                    %AA=Vel.depth_xducer;
                    Vel.depth_xducer(ibb)=3276.8.*2+Vel.depth_xducer(ibb);
                end
                
                
                %Set GUI's
                prdown_ff=1.887;%4933/2617; %empirical fudge factor for adjusting ADCP pressure based on CTD pressure
                prdown_ff=1; %downlooker now has correct cal so correction no longer needed.
                set(handles.ensemble_text,'String',num2str(Vel.ens_no));
                set(handles.pressure_text,'String',num2str(0.1*round(Vel.depth_xducer*prdown_ff*10))); %plot only decimal place
                set(handles.pitch_text,'String',num2str(Vel.pitch));
                set(handles.roll_text,'String',num2str(Vel.roll));
                set(handles.heading_text,'String',num2str(Vel.heading));
                %            set(handles.range_text,'String',num2str(.1*round(Vel.btrange(1)/100*10))); %put in m; round to 1 sig fig
                set(handles.range_text,'String',num2str(.1*round(nanmean(Vel.btrange)/100*10))); %put in m; round to 1 sig fig; take nanmin 1/2014
                set(handles.datestr_text,'String',datestr(Vel.dtnum));
                % set(handles.batt_text,'String',num2str(Vel.adcchan1*0.27));
                set(handles.batt_text,'String',num2str(Vel.adcchan1*0.36176));
            else
                set(handles.ensemble_text,'String','-');
                set(handles.pressure_text,'String','-');
                set(handles.pitch_text,'String','-');
                set(handles.roll_text,'String','-');
                set(handles.heading_text,'String','-');
                set(handles.range_text,'String','-');
                set(handles.datestr_text,'String','-')
                set(handles.batt_text,'String','-')
            end
            
            %plot
            %    plot(handles.plot1_axes,Vel.ec1_bm,Vel.z_adcp,...
            %        handles.plot1_axes,Vel.ec2_bm,Vel.z_adcp,...
            %        handles.plot1_axes,Vel.ec3_bm,Vel.z_adcp,...
            %        handles.plot1_axes,Vel.ec4_bm,Vel.z_adcp)
            if ~isempty(Vel)
                plot(handles.plot1_axes,[Vel.ec1_bm Vel.ec2_bm Vel.ec3_bm Vel.ec4_bm],Vel.z_adcp)
            else
                plot(handles.plot1_axes,nan,nan)
            end
            set(handles.plot1_axes,'ydir','reverse')
            set(handles.plot1_axes,'xlim',[0 256])
            set(handles.plot1_axes,'xgrid','on','ygrid','on')
            set(handles.plot1_axes,'ylim',[0 200])
            ylabel(handles.plot1_axes,'z / m')
            xlabel(handles.plot1_axes,'Echo Intensity')
            
            
            
            %1/28/2014:
            %Check if we are in range of the bottom.
            if handles.enable_bottom_warnings
                if ~isempty(Vel)
                    
                    %%
                    thisrange=nanmean(Vel.btrange)/100;
                    if thisrange < 120 % ~isnan(thisrange)  %if last ping we saw the bottom and <100m range
                        if thisrange < handles.btrange
                            handles.decreasing=1;
                        else
                            handles.decreasing=0;
                        end
                        
                        handles.btrange= thisrange; %set the new bottom range - average of all four beams, in m
                        handles.bthits=handles.bthits+1;
                        
                    else %if a nan then reset
                        handles.decreasing=0;
                        handles.bthits=0;
                        
                    end
                    
                    if handles.bthits >=handles.hits_to_warn & handles.decreasing & now > handles.last + handles.max_warn_freq_sec/3600/24
                        %get people's attention...
                        %%
                        WaveCell={'mistake_hg.wav','vibe.wav','pay.wav','whatchudoin.wav','garycoleman.wav','dontgive'};
                        mrt=imread('MrT.jpg');
                        ri=randi(length(WaveCell));
                        %ri=1;
                        Y4=wavread([WaveCell{ri}]);
                        %            Y4=wavread(['inhole1.wav']);
                        %            Y4=wavread(['2turds.wav']);
                        P4=audioplayer(Y4,11025);
                        f = figure;
                        image(mrt);set(gca,'xtick',[],'ytick',[])
                        play(P4)
                        
                        title('Pay attention! Bottom approaching!')
                        h = uicontrol('Position', [20 20 200 40], 'String', 'I''m awake', ...
                            'Callback', 'uiresume(gcbf)');
                        uiwait(gcf);
                        handles.last=now;
                        close(f);
                    end
                end
            end %end if bottom warnings enabled
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Now new code for v1.1 - add to a buffer of Vel and plot in a new window
            if handles.enable_plots %enable new code
                if ~isempty(Vel) %if we have valid data
                    
                    %fields to copy directly from the Vel structure
                    fields={'depth_xducer';'yday';'dtnum';'heading';'pitch';'roll';'btrange';'ec1_bm';'u';'v';'btvel_bm';'ens_no'};
                    %fields we need to interpolate onto our depth grid
                    fields2={'vb';'vb_bt';'ub';'ub_bt'};
                    
                    if ~isfield(handles,'Vela') %make Vela if we have not yet created it.
                        %Create Vela with the following fields, straight from Vel
                        %Do this the first time we have a valid Vel - this way we know the correct
                        %size.
                        handles.Vela.z_adcp=Vel.z_adcp;
                        for c=1:length(fields)
                            [m,n]=size(Vel.(fields{c}));
                            handles.Vela.(fields{c})=nan(m,handles.no_ensembles);
                        end
                        
                        %Then open the plotting window
                        figure(1)
                        set(gcf,'position',[740    46   536   677])
                        %figure(fig) %always get back to the GUI window after plotting in other windows
                    end
                    
                    if handles.start_ensemble==1 %if we need to create grids (ie if plotting range changed).
                        %Now make the depth vector
                        handles.Vela.zg=(handles.zmin:4:handles.zmax)';
                        %Then make the fields.  This gets called when we start, or when we change
                        %the grid
                        for c=1:length(fields2)
                            handles.Vela.(fields2{c})=nan(length(handles.Vela.zg),handles.no_ensembles);
                        end
                    end
                    
                    %Now add the needed fields of Vel to Vela, removing
                    %the first record.
                    for c=1:length(fields)
                        handles.Vela.(fields{c})=[handles.Vela.(fields{c})(:,2:end) Vel.(fields{c})];
                    end
                    
                    %shift the gridded fields over by one - then we'll add the last
                    %one in.
                    for c=1:length(fields2)
                        [m,n]=size(handles.Vela.(fields2{c}));
                        handles.Vela.(fields2{c})=[handles.Vela.(fields2{c})(:,2:end) nan(m,1)];
                    end
                    
                    %Interpolate the new data onto the grid
                    ido=handles.start_ensemble:handles.no_ensembles;
                    ig=find(~isnan(handles.Vela.yday(ido)));
                    for c=1:length(ig)
                        igz=find(~isnan(handles.Vela.u(:,ido(ig(c)))));
                        if length(igz) >=2
                            handles.Vela.ub(:,ido(ig(c)))=interp1(prdown_ff*handles.Vela.depth_xducer(ido(ig(c))) + handles.Vela.z_adcp(igz),handles.Vela.u(igz,ido(ig(c))),handles.Vela.zg);
                            handles.Vela.vb(:,ido(ig(c)))=interp1(prdown_ff*handles.Vela.depth_xducer(ido(ig(c))) + handles.Vela.z_adcp(igz),handles.Vela.v(igz,ido(ig(c))),handles.Vela.zg);
                            
                            handles.Vela.ub_bt(:,ido(ig(c)))=interp1(prdown_ff*handles.Vela.depth_xducer(ido(ig(c))) + handles.Vela.z_adcp(igz),...
                                handles.Vela.u(igz,ido(ig(c)))-handles.Vela.btvel_bm(1,ido(ig(c))),handles.Vela.zg);
                            handles.Vela.vb_bt(:,ido(ig(c)))=interp1(prdown_ff*handles.Vela.depth_xducer(ido(ig(c))) + handles.Vela.z_adcp(igz),...
                                handles.Vela.v(igz,ido(ig(c)))-handles.Vela.btvel_bm(2,ido(ig(c))),handles.Vela.zg);
                        end
                        %igz=find(handles.Vela.ec1_bm(:,ido(ig(c))));
                        %if length(igz) >=2
                        %    handles.Vela.ec1b(:,ido(ig(c)))=interp1(prdown_ff*handles.Vela.depth_xducer(ido(ig(c))) + handles.Vela.z_adcp(igz),...
                        %        handles.Vela.ec1_bm(igz,ido(ig(c))),handles.Vela.zg);
                        %end
                    end
                    
                    %indicate that we only need to do the next one
                    handles.start_ensemble=handles.no_ensembles;
                    
                    %update plotting range for the next one if we are out of it
                    if Vel.depth_xducer*prdown_ff < handles.zmin
                        old_diff=handles.zmax - handles.zmin;
                        handles.zmax = fix(Vel.depth_xducer*prdown_ff) + old_diff/2;
                        handles.zmin = handles.zmax -old_diff;
                        set(handles.edit_zmin,'String',num2str(handles.zmin));
                        set(handles.edit_zmax,'String',num2str(handles.zmax));
                        handles.start_ensemble=1;
                    end
                    if Vel.depth_xducer*prdown_ff > handles.zmax
                        old_diff=handles.zmax - handles.zmin;
                        handles.zmin = fix(Vel.depth_xducer*prdown_ff) - old_diff / 2;
                        handles.zmax = handles.zmin +old_diff;
                        set(handles.edit_zmin,'String',num2str(handles.zmin));
                        set(handles.edit_zmax,'String',num2str(handles.zmax));
                        handles.start_ensemble=1;
                        
                    end
                    
                    %                plot_recent_data(handles)
                    ig=find(~isnan(handles.Vela.dtnum));
                    if length(ig) >=3 & rem(handles.count,5)==0 %only plot every so often
                        use_ezcf=1;
                        sub_imagesc=1;
                        use_image=0;
                        do_shear=1;
                        
                        dz = diff(handles.Vela.zg(1:2));
                        nbins = size(handles.Vela.ub,1);
                        nprofs = length(ig);
                        
                        umeas = handles.Vela.ub(:,ig);
                        vmeas = handles.Vela.vb(:,ig);
                        
                        % Mask below bottom return
                        zbot = handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100;
                        bd = find(repmat(handles.Vela.zg,1,nprofs)>repmat(zbot,nbins,1));
                        umeas(bd) = NaN;
                        vmeas(bd) = NaN;
                        
                        upack = repmat(0,1,nprofs);
                        vpack = repmat(0,1,nprofs);
                        
                        if do_shear & nprofs>10,
                            % compute shear solution
                            
                            for ii=1:size(umeas,1),
                                % filter in time
                                %umeas(ii,:) = medfilt(umeas(ii,:),5);
                                %vmeas(ii,:) = medfilt(vmeas(ii,:),5);
                                
                                umeas(ii,:) = stdfilt(umeas(ii,:),9,1);
                                vmeas(ii,:) = stdfilt(vmeas(ii,:),9,1);
                                
                            end
                            
                            dudz = -diff(umeas)/dz;
                            dvdz = -diff(vmeas)/dz;
                            
                            Su = nanmean(dudz')';
                            Sv = nanmean(dvdz')';
                            ussol = repmat(NaN,nbins,1);
                            vssol = repmat(NaN,nbins,1);
                            gd = find(isfinite(Su)&isfinite(Sv));
                            ussol(gd) = -cumsum(Su(gd))*dz; % single u profile
                            vssol(gd) = -cumsum(Sv(gd))*dz;
                            % constrain mean profile over grid to be zero.
                            % alternatively could set a particular depth to be
                            % zero.
                            ussol = ussol-nanmean(ussol);
                            vssol = vssol-nanmean(vssol);
                            
                            %size(Su)
                            %size(ussol)
                            %size(handles.Vela.ub)
                            %size(ig)
                            
                            % Now match individual profiles with shear solution
                            % in overlap range to get package motion (relative to
                            % zero point)
                            %upack = nanmean(repmat(ussol,1,nprofs) - handles.Vela.ub(:,ig));
                            %vpack = nanmean(repmat(vssol,1,nprofs) - handles.Vela.vb(:,ig));
                            upack = nanmean(repmat(ussol,1,nprofs) - umeas);
                            vpack = nanmean(repmat(vssol,1,nprofs) - vmeas);
                            
                            
                            % Check whether valid bottom track data exists and (a)
                            % shift upack,vpack to agree over time of overlap and
                            % (b) substitute bottom track vel in overlap region
                            udiff = upack + handles.Vela.btvel_bm(1,ig);
                            vdiff = vpack + handles.Vela.btvel_bm(2,ig);
                            %mudiff = nanmean(udiff);
                            %mvdiff = nanmean(vdiff);
                            mudiff = robmean(udiff,0.2);
                            mvdiff = robmean(vdiff,0.2)
                            max(vdiff)
                            min(vdiff)
                            
                            isbt = find(isfinite(handles.Vela.btvel_bm(1,ig)));
                            is3 = find(diff(diff(isbt))==0); % 3 consecutive BT
                            
                            %if isfinite(mudiff) & isfinite(mvdiff),
                            if ~isempty(is3)&isfinite(mudiff)&isfinite(mvdiff),
                                % enough bottom tracking exists, so apply corrections
                                upack = upack - mudiff;
                                vpack = vpack - mvdiff;
                                
                                %upack(isbt) = -handles.Vela.btvel_bm(1,ig(isbt));
                                %vpack(isbt) = -handles.Vela.btvel_bm(2,ig(isbt));
                            end
                            
                            % and filter package velocity before adding back in
                            upack = stdfilt(upack,9,1);
                            vpack = stdfilt(vpack,9,1);
                            
                        end % if do_shear & nprofs>10
                        
                        if use_ezcf
                            % Plot - a dtnum version
                            figure(1)
                            clf
                            cl=[handles.vmin handles.vmax];
                            cl_i=[0 150];
                            yl=[handles.zmin handles.zmax];
                            ncolors=16;
                            %                xl=[min(handles.Vela.dtnum(ig)) max(handles.Vela.dtnum(ig))];
                            dt=handles.no_ensembles * 1.5 * 1.4 / 3600 / 24; %compute the window aperture from #ensembles and ping rate
                            xl=max(handles.Vela.dtnum(ig))+[-dt 0];
                            ax1=MySubplot(.1,0.04,0,0.05,.67,0.01,1,3);
                            ax2=MySubplot(.1,0.04,0.02,.37,.07,0.02,2,2);
                            %
                            axes(ax1(1));
                            plot(handles.Vela.dtnum(ig),handles.Vela.heading(ig));
                            xlim(xl)
                            datetick('x','keeplimits')
                            xtloff
                            SubplotLetter('Heading',.01,.8)
                            ylim([0 360])
                            ylabel('degrees')
                            set(gca,'ytick',[0 90 180 270 360])
                            grid
                            
                            axes(ax1(2))
                            plot(handles.Vela.dtnum(ig),handles.Vela.pitch(ig), 'b-',handles.Vela.dtnum(ig),handles.Vela.roll(ig),'g-');
                            xlim(xl)
                            datetick('x','keeplimits')
                            xtloff
                            SubplotLetter('pitch (blue), roll (green)',.01,.9)
                            grid
                            ylabel('degrees')
                            axes(ax1(3))
                            %
                            if sub_imagesc
                                imagesc(handles.Vela.dtnum(ig),handles.Vela.z_adcp,handles.Vela.ec1_bm(:,ig),cl_i);
                            else
                                ezcf(handles.Vela.dtnum(ig),handles.Vela.z_adcp,handles.Vela.ec1_bm(:,ig),cl_i,ncolors);
                            end
                            %
                            %
                            %caxis([-.25 .25])
                            xlim(xl)
                            datetick('x','keeplimits')
                            %xtloff
                            %
                            grid on
                            ylabel('Range / m')
                            SubplotLetter('Beam 1 intensity',.01,.88)
                            axes(ax2(1))
                            %ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.ub(:,ig),cl,ncolors)
                            %ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.ub(:,ig)+repmat(upack,nbins,1),cl,ncolors)
                            if sub_imagesc
                                pdata = umeas+repmat(upack,nbins,1);
                                bd = find(isnan(pdata));
                                pdata(bd) = 0;
                                imagesc(handles.Vela.dtnum(ig),handles.Vela.zg,pdata,cl)
                            else
                                ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,umeas+repmat(upack,nbins,1),cl,ncolors)
                            end
                            colormap(redblue2(32))
                            
                            hold on
                            plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-')
                            plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                            xlim(xl)
                            grid on
                            datetick('x','keeplimits')
                            xtloff
                            ylim(yl)
                            SubplotLetter('Measured u',.01,.9)
                            ylabel('P / dbar')
                            axes(ax2(2))
                            
                            %ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.vb(:,ig),cl,ncolors)
                            %ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.vb(:,ig)+repmat(vpack,nbins,1),cl,ncolors)
                            if sub_imagesc
                                pdata = vmeas+repmat(vpack,nbins,1);
                                bd = find(isnan(pdata));
                                pdata(bd) = 0;
                                imagesc(handles.Vela.dtnum(ig),handles.Vela.zg,pdata,cl)
                            else
                                ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,vmeas+repmat(vpack,nbins,1),cl,ncolors)
                            end
                            hold on;
                            plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-')
                            plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                            xlim(xl)
                            ylim(yl)
                            grid on
                            datetick('x','keeplimits')
                            SubplotLetter('Measured v',.01,.9)
                            %ylabel('P / dbar')
                            xtloff
                            xyloff
                            
                            axes(ax2(3))
                            if length(find(~isnan(nanmean(handles.Vela.ub_bt(:,ig))))) > 5
                                if sub_imagesc
                                    pdata = handles.Vela.ub_bt(:,ig);
                                    bd = find(isnan(pdata));
                                    pdata(bd) = 0;
                                    imagesc(handles.Vela.dtnum(ig),handles.Vela.zg,pdata,cl)
                                else
                                    ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.ub_bt(:,ig),cl,ncolors)
                                end
                            end
                            hold on
                            plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-')
                            plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                            xlim(xl)
                            ylim(yl)
                            grid on
                            datetick('x','keeplimits')
                            %xtloff
                            SubplotLetter('u w/ bottom track',.01,.9)
                            xlabel('Time UTC')
                            ylabel('P / dbar')
                            axis ij
                            axes(ax2(4))
                            
                            if length(find(~isnan(nanmean(handles.Vela.ub_bt(:,ig))))) > 5
                                if sub_imagesc
                                    pdata = handles.Vela.vb_bt(:,ig);
                                    bd = find(isnan(pdata));
                                    pdata(bd) = 0;
                                    imagesc(handles.Vela.dtnum(ig),handles.Vela.zg,pdata,cl)
                                else
                                    ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.vb_bt(:,ig),cl,ncolors)
                                end
                            end
                            hold on
                            plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-')
                            plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                            xlim(xl)
                            ylim(yl)
                            grid on
                            datetick('x','keeplimits')
                            SubplotLetter('v w/ bottom track',.01,.9)
                            %ylabel('P / dbar')
                            xyloff
                            axis ij
                            xlabel('Time UTC')
                            %colorbar
                        elseif use_image
                            % attempt to get datenum without redrawing
                            %disp(['start_ens = ' num2str(handles.start_ensemble)]);
                            do_datetick = 0;
                            figure(1)
                            if handles.start_ensemble==1,
                                % First ensemble: draw all plots from scratch
                                % and put handles in p structure
                                
                                clf
                                
                                cl=[handles.vmin handles.vmax];
                                cl_i=[0 150];
                                yl=[handles.zmin handles.zmax];
                                ncolors=16;
                                %                xl=[min(handles.Vela.dtnum(ig)) max(handles.Vela.dtnum(ig))];
                                dt=handles.no_ensembles * 1.5 * 1.4 / 3600 / 24; %compute the window aperture from #ensembles and ping rate
                                xl=max(handles.Vela.dtnum(ig))+[-dt 0];
                                ax=MySubplot(.1,.1,0,.1,.1,0.02,1,7);
                                %
                                % First plot: heading
                                axes(ax(1));
                                p.h1=plot(handles.Vela.dtnum(ig),handles.Vela.heading(ig));
                                xlim(xl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                xtloff
                                SubplotLetter('Heading',.01,.9)
                                ylim([0 360])
                                ylabel('degrees')
                                set(gca,'ytick',[0 90 180 270 360])
                                grid
                                
                                % second plot: pitch and roll
                                axes(ax(2))
                                p.h2=plot(handles.Vela.dtnum(ig),handles.Vela.pitch(ig), 'b-',handles.Vela.dtnum(ig),handles.Vela.roll(ig),'g-');
                                xlim(xl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                xtloff
                                SubplotLetter('pitch (blue), roll (green)',.01,.9)
                                grid
                                ylabel('degrees')
                                
                                % third plot: beam 1 echo strength
                                axes(ax(3))
                                %
                                %ezcf(handles.Vela.dtnum(ig),handles.Vela.z_adcp,handles.Vela.ec1_bm(:,ig),cl_i,ncolors);
                                p.h3=imagesc(handles.Vela.dtnum(ig),handles.Vela.z_adcp,handles.Vela.w(:,ig),cl_i);
                                %
                                %
                                %caxis([-.25 .25])
                                xlim(xl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                xtloff
                                %
                                ylabel('Range / m')
                                SubplotLetter('Beam 1 intensity',.01,.9)
                                
                                % fourth plot: u relative
                                axes(ax(4))
                                %ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.ub(:,ig),cl,ncolors)
                                p.h4=imagesc(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.ub(:,ig),cl)
                                hold on
                                p.h4a=plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-');
                                p.h4b=plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks');
                                end
                                xtloff
                                ylim(yl)
                                SubplotLetter('Measured u',.01,.9)
                                ylabel('P / dbar')
                                
                                % fifth plot: v relative
                                axes(ax(5))
                                
                                %ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.vb(:,ig),cl,ncolors)
                                %size(handles.Vela.dtnum(ig))
                                %size(handles.Vela.zg)
                                %size(handles.Vela.vb(:,ig))
                                %Note: colormap has 64 values, but ncolors is 16
                                p.h5=imagesc(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.vb(:,ig),cl)
                                hold on; caxis(cl);
                                p.h5a=plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-');
                                p.h5b=plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl);
                                ylim(yl);
                                if do_datetick
                                    datetick('x','keeplimits','keepticks');
                                end
                                
                                SubplotLetter('Measured v',.01,.9)
                                ylabel('P / dbar')
                                xtloff
                                
                                % sixth plot: bottom-tracked u
                                axes(ax(6))
                                if length(find(~isnan(nanmean(handles.Vela.ub_bt(:,ig))))) > 5
                                    %ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.ub_bt(:,ig),cl,ncolors)
                                    p.h6=imagesc(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.ub_bt(:,ig),cl);
                                end
                                hold on
                                p.h6a=plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-')
                                p.h6b=plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl)
                                ylim(yl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                xtloff
                                SubplotLetter('u w/ bottom track',.01,.9)
                                ylabel('P / dbar')
                                axis ij
                                
                                % seventh plot: bottom-tracked v
                                axes(ax(7))
                                
                                if length(find(~isnan(nanmean(handles.Vela.ub_bt(:,ig))))) > 5
                                    %ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.vb_bt(:,ig),cl,ncolors)
                                    p.h7=imagesc(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.vb_bt(:,ig),cl);
                                end
                                hold on
                                p.h7a=plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-');
                                p.h7b=plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl)
                                ylim(yl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                SubplotLetter('v w/ bottom track',.01,.9)
                                ylabel('P / dbar')
                                axis ij
                                xlabel('Time UTC')
                                %colorbar
                                
                            else % if not first ensemble, just update plots
                                
                                cl=[handles.vmin handles.vmax];
                                cl_i=[0 150];
                                yl=[handles.zmin handles.zmax];
                                ncolors=16;
                                %                xl=[min(handles.Vela.dtnum(ig)) max(handles.Vela.dtnum(ig))];
                                dt=handles.no_ensembles * 1.5 * 1.4 / 3600 / 24; %compute the window aperture from #ensembles and ping rate
                                xl=max(handles.Vela.dtnum(ig))+[-dt 0];
                                ax=MySubplot(.1,.1,0,.1,.1,0.02,1,7);
                                %
                                axes(ax(1));
                                p.h1=plot(handles.Vela.dtnum(ig),handles.Vela.heading(ig));
                                xlim(xl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                xtloff
                                SubplotLetter('Heading',.01,.9)
                                ylim([0 360])
                                ylabel('degrees')
                                set(gca,'ytick',[0 90 180 270 360])
                                grid
                                axes(ax(2))
                                plot(handles.Vela.dtnum(ig),handles.Vela.pitch(ig), 'b-',handles.Vela.dtnum(ig),handles.Vela.roll(ig),'g-');
                                xlim(xl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                xtloff
                                SubplotLetter('pitch (blue), roll (green)',.01,.9)
                                grid
                                ylabel('degrees')
                                axes(ax(3))
                                %
                                ezcf(handles.Vela.dtnum(ig),handles.Vela.z_adcp,handles.Vela.ec1_bm(:,ig),cl_i,ncolors);
                                %
                                %
                                %caxis([-.25 .25])
                                xlim(xl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                xtloff
                                %
                                ylabel('Range / m')
                                SubplotLetter('Beam 1 intensity',.01,.9)
                                axes(ax(4))
                                ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.ub(:,ig),cl,ncolors)
                                %image(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.ub(:,ig),cl,ncolors)
                                hold on
                                plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-')
                                plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                xtloff
                                ylim(yl)
                                SubplotLetter('Measured u',.01,.9)
                                ylabel('P / dbar')
                                axes(ax(5))
                                
                                %ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.vb(:,ig),cl,ncolors)
                                %size(handles.Vela.dtnum(ig))
                                %size(handles.Vela.zg)
                                %size(handles.Vela.vb(:,ig))
                                %Note: colormap has 64 values, but ncolors is 16
                                imagesc(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.vb(:,ig),cl)
                                hold on; caxis(cl);
                                plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-')
                                plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl)
                                ylim(yl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                
                                SubplotLetter('Measured v',.01,.9)
                                ylabel('P / dbar')
                                xtloff
                                
                                axes(ax(6))
                                if length(find(~isnan(nanmean(handles.Vela.ub_bt(:,ig))))) > 5
                                    ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.ub_bt(:,ig),cl,ncolors)
                                end
                                hold on
                                plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-')
                                plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl)
                                ylim(yl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                xtloff
                                SubplotLetter('u w/ bottom track',.01,.9)
                                ylabel('P / dbar')
                                axis ij
                                axes(ax(7))
                                
                                if length(find(~isnan(nanmean(handles.Vela.ub_bt(:,ig))))) > 5
                                    ezcf(handles.Vela.dtnum(ig),handles.Vela.zg,handles.Vela.vb_bt(:,ig),cl,ncolors)
                                end
                                hold on
                                plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff,'k-')
                                plot(handles.Vela.dtnum(ig),handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl)
                                ylim(yl)
                                if do_datetick
                                    datetick('x','keeplimits','keepticks')
                                end
                                SubplotLetter('v w/ bottom track',.01,.9)
                                ylabel('P / dbar')
                                axis ij
                                xlabel('Time UTC')
                                %colorbar
                                
                            end
                            
                            
                        else %use pcolor, but it doesn't work with dtnum
                            time=(handles.Vela.yday(ig) - handles.Vela.yday(ig(end)))/24/60;
                            % Plot - a dtnum version
                            ig=find(~isnan(handles.Vela.dtnum));
                            if length(ig) >=4
                                figure(1)
                                clf
                                cl=[-.45 .45];
                                cl_i=[0 200];
                                ncolors=16;
                                %                xl=[min(handles.Vela.dtnum(ig)) max(handles.Vela.dtnum(ig))];
                                dt=handles.no_ensembles * 1.4 / 3600 / 24; %compute the window aperture from #ensembles and ping rate
                                xl=[-dt*60 0];
                                ax=MySubplot(.1,.1,0,.1,.1,0.02,1,7)
                                %
                                axes(ax(1))
                                plot(time,handles.Vela.heading(ig));
                                xlim(xl)
                                %datetick('x','keeplimits','keepticks')
                                xtloff
                                SubplotLetter('Heading',.01,.9)
                                ylim([0 360])
                                set(gca,'ytick',[0 90 180 270 360])
                                grid
                                axes(ax(2))
                                plot(time,handles.Vela.pitch(ig), 'b-',time,handles.Vela.roll(ig),'g-');
                                xlim(xl)
                                %datetick('x','keeplimits','keepticks')
                                xtloff
                                SubplotLetter('pitch (blue), roll (green)',.01,.9)
                                grid
                                axes(ax(3))
                                %
                                ezpc(time,handles.Vela.z_adcp,handles.Vela.ec1_bm(:,ig));
                                %
                                hold on
                                plot(time,handles.Vela.depth_xducer(ig)*prdown_ff);
                                hold off
                                %
                                %caxis([-.25 .25])
                                xlim(xl)
                                %datetick('x','keeplimits','keepticks')
                                xtloff
                                %
                                ylabel('Range / m')
                                SubplotLetter('Beam 1 intensity',.01,.9)
                                %
                                axes(ax(4))
                                ezpc(time,handles.Vela.zg,handles.Vela.ub(:,ig))
                                hold on
                                plot(time,handles.Vela.depth_xducer(ig)*prdown_ff,'k-')
                                plot(time,handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl)
                                %datetick('x','keeplimits','keepticks')
                                xtloff
                                SubplotLetter('Measured u',.01,.9)
                                ylabel('P / dbar')
                                axes(ax(5))
                                
                                ezpc(time,handles.Vela.zg,handles.Vela.vb(:,ig))
                                hold on
                                plot(time,handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl)
                                %datetick('x','keeplimits','keepticks')
                                SubplotLetter('Measured v',.01,.9)
                                ylabel('P / dbar')
                                xtloff
                                
                                axes(ax(6))
                                ezpc(time,handles.Vela.zg,handles.Vela.ub_bt(:,ig))
                                hold on
                                plot(time,handles.Vela.depth_xducer(ig)*prdown_ff,'k-')
                                plot(time,handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl)
                                %datetick('x','keeplimits','keepticks')
                                xtloff
                                SubplotLetter('u w/ bottom track',.01,.9)
                                ylabel('P / dbar')
                                axes(ax(7))
                                
                                ezpc(time,handles.Vela.zg,handles.Vela.vb_bt(:,ig))
                                hold on
                                plot(time,handles.Vela.depth_xducer(ig)*prdown_ff+nanmean(handles.Vela.btrange(:,ig))/100,'k-');
                                xlim(xl)
                                %                    datetick('x','keeplimits','keepticks')
                                SubplotLetter('v w/ bottom track',.01,.9)
                                ylabel('P / dbar')
                                
                                xlabel('Time (min)')
                            end
                        end
                        %                    figure(fig) %always get back to the GUI window after plotting in other windows
                    end %if length(ig) > 10
                end %end if exist(vel)
            end %if 0
        end %end if rem(count,skip)
    end %if no mod errors
end %end if complete record


% Update handles structure
guidata(fig, handles);


% --- Executes on button press in stop_button.
function stop_button_Callback(hObject, eventdata, handles)
% hObject    handle to stop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%err=fclose(handles.serial_fid);
fclose(handles.serial_s);

figure(1)
close(gcf)

delete('tmpfil.000')

%return status to idle.
handles.status='Idle.';
set(handles.status_text,'String',handles.status)

%set(handles.filenum,'String','');

%close the file if open.
if handles.out >0
    err=fclose(handles.out);
end

%turn the outfilename and directory edit boxes back on
set(handles.outfilename_edit,'Enable','On');
set(handles.outputdir_edit,'Enable','On');

%set button enable status
set(handles.start_button,'Enable','On');
set(handles.stop_button,'Enable','Off');

% Update handles structure
guidata(hObject, handles);

function outputdir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to outputdir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputdir_edit as text
%        str2double(get(hObject,'String')) returns contents of outputdir_edit as a double

handles.outputdir=get(hObject,'String');

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function outputdir_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outputdir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Outputs from this function are returned to the command line.
function varargout = LADCP_DAQ_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
figure(1)
close(gcf)
varargout{1} = handles.output;



function outfilename_edit_Callback(hObject, eventdata, handles)
% hObject    handle to outfilename_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outfilename_edit as text
%        str2double(get(hObject,'String')) returns contents of outfilename_edit as a double

handles.outfilename=get(hObject,'String');
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function outfilename_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outfilename_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_zmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zmax as text
%        str2double(get(hObject,'String')) returns contents of edit_zmax as a double


%handles.zmax=str2num(get(hObject,'String'));
%handles.start_ensemble=1;

%Don't allow setting less than zmin
val=str2num(get(hObject,'String'));

if val <= handles.zmin
    set(hObject,'String',num2str(handles.zmax))
else
    handles.zmax=val;
    handles.start_ensemble=1;
end

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_zmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.zmax=str2num(get(hObject,'String'));
handles.start_ensemble=1;
guidata(hObject,handles)



function edit_zmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zmin as text
%        str2double(get(hObject,'String')) returns contents of edit_zmin as a double

%Don't allow setting greater than zmax
val=str2num(get(hObject,'String'));

if val >= handles.zmax
    set(hObject,'String',num2str(handles.zmin))
else
    handles.zmin=val;
    handles.start_ensemble=1;
    
end

%handles.zmin=str2num(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_zmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.zmin=str2num(get(hObject,'String'));
handles.start_ensemble=1;
guidata(hObject,handles)



function edit_vmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vmax as text
%        str2double(get(hObject,'String')) returns contents of edit_vmax as a double

%Don't allow setting less than vmin
val=str2num(get(hObject,'String'));

if val <= handles.vmin
    set(hObject,'String',num2str(handles.vmax))
else
    handles.vmax=val;
    
end

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_vmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_vmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vmin as text
%        str2double(get(hObject,'String')) returns contents of edit_vmin as a double

%handles.vmin=str2num(get(hObject,'String'));

%Don't allow setting greater than vmax
val=str2num(get(hObject,'String'));

if val >= handles.vmax
    set(hObject,'String',num2str(handles.vmin))
else
    handles.vmin=val;
    
end

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_vmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_plot.
function checkbox_plot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot

handles.enable_plots=get(hObject,'Value');

guidata(hObject,handles)
