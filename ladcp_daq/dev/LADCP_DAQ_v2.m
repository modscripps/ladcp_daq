%LADCP_DAQ.m
%Version 1.0, MHA
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

% Last Modified by GUIDE v2.5 25-Jul-2012 02:05:49

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
handles.outputdir='C:\Documents and Settings\All Users\Desktop\LADCP_DAQ\data\';
handles.outfilename='RR1209_LADCP_Downlooker_cast1';

%initialize counters etc
handles.bufferbytes=0; %bytes in serial buffer
handles.filebytes=0; %bytes written to file
handles.lastline='';
handles.status='Idle.';

handles.filenum=0; %increment this after a certain number of records
handles.recsperfile=100; %number of records to store per file
handles.MAXFILEBYTES=5e6; %size of each file in bytes. 5 MB
handles.plotskip=3; %plot every (this)'th record; can set to greater than 1 if processing can't keep up

%set GUI widgets to initial values
set(handles.outfilename_edit,'String',handles.outfilename);
set(handles.outputdir_edit,'String',handles.outputdir);
set(handles.lastline_text,'String',handles.lastline);
set(handles.status_text,'String',handles.status);

set(handles.bufferbytes_text,'String',num2str(handles.bufferbytes));
set(handles.filebytes_text,'String',num2str(handles.filebytes));

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
i1=findstr(handles.buffa,'7F7F93');

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

    %Then turn this into an unsigned integer by calling hex2dec on each column
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
        %Original version: seems odd with pressure ~ real pres / 2
        Vel=Get_ADCP_fullSC_BT('tmpfil.000',-1,1);
        %John's "new" version for LADCP
%        Vel=Get_ADCP_fullSC_LADCP('tmpfil.000',-1,1);
        if ~isempty(Vel)
            %Set GUI's
            fudgefactor=1.887;%4933/2617; %empirical fudge factor for adjusting ADCP pressure based on CTD pressure
            set(handles.ensemble_text,'String',num2str(Vel.ens_no));
            set(handles.pressure_text,'String',num2str(Vel.depth_xducer*fudgefactor));
            set(handles.pitch_text,'String',num2str(Vel.pitch));
            set(handles.roll_text,'String',num2str(Vel.roll));
            set(handles.heading_text,'String',num2str(Vel.heading));
            set(handles.range_text,'String',num2str(nanmin(Vel.btrange)/100)); %put in m
            %set(handles.range_text,'String',num2str(Vel.bottomBT/100));
             %%put in m
            set(handles.datestr_text,'String',datestr(Vel.dtnum));
        else
            set(handles.ensemble_text,'String','-');
            set(handles.pressure_text,'String','-');
            set(handles.pitch_text,'String','-');
            set(handles.roll_text,'String','-');
            set(handles.heading_text,'String','-');
            set(handles.range_text,'String','-');
            set(handles.datestr_text,'String','-')
        end

        %plot
        %    plot(handles.plot1_axes,Vel.ec1_bm,Vel.z_adcp,...
        %        handles.plot1_axes,Vel.ec2_bm,Vel.z_adcp,...
        %        handles.plot1_axes,Vel.ec3_bm,Vel.z_adcp,...
        %        handles.plot1_axes,Vel.ec4_bm,Vel.z_adcp)
        if ~isempty(Vel)
            plot(handles.plot1_axes,Vel.ec1_bm,Vel.z_adcp)
        else
            plot(handles.plot1_axes,nan,nan)
        end
        
        if ~isempty(Vel)
            if ~isnan(Vel.btrange(1));
            plot(handles.plot2_axes,Vel.u+Vel.btvel_bm(1,:),Vel.z_adcp+Vel.depth_xducer,'b'); %can't recall if plus or minus..
            hold on
            plot(handles.plot2_axes,Vel.v+Vel.btvel_bm(2,:),Vel.z_adcp+Vel.depth_xducer,'r');
            hold off
            else
            plot(handles.plot2_axes,Vel.v,Vel.z_adcp+Vel.depth_xducer,'r')
            hold on
            plot(handles.plot2_axes,Vel.u,Vel.z_adcp+Vel.depth_xducer,'b')
            %legend(handles.plot2_axes,'v (red)','u (blue)');
            hold off
            end
             
        else
            plot(handles.plot2_axes,nan,nan);
        end
        
        
            
            
        set(handles.plot1_axes,'ydir','reverse')
        set(handles.plot1_axes,'xlim',[0 256])
        set(handles.plot1_axes,'ylim',[0 400])
        ylabel(handles.plot1_axes,'z / m')
        Xlabel(handles.plot1_axes,'Echo Intensity')
    end

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

