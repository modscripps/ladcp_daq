function raw2mat_SC(transname,matbasename,inrule,outrule,inrange,ens_ind_rng);
%
% RAW2MAT_SC translates raw SC files and outputs a series of matlab 
% 
% USAGE: raw2mat_SC(transname,matbasename,inrule,outrule,inrange, ens_ind_rng);
% 
%  This program translates a set of transect files.  It only requires 
%  the prefix from the raw files to run; it incrementally will
%  read all the files it finds in the same directory.  i.e. if you
%  have a deployment MC97001.000, MC97001.001, MC97001.002,...
%  then a call like raw2mat('MC97001',MC97001_') will create
%  a series of files MC97001_1.mat MC97001_2.mat, ....  The numbering
%  of the output files does not correspond to the numbering of the
%  TRANSECT files.  A new matlab file is written every 900 ensembles.

% ens_ind_rng = min,max range for indices of ensembles in file (1 to ???)
%
%  Almost all the data supplied by the ADCP is translated to 
%  matlab variables.  A couple of things have been left out.
%

%  Copyright: Jody Klymak, Nov 19. 1998: jklymak@apl.washington.edu
%  Conditions:  This routine is supplied with no warranty.  It is freely
%  distributable under the condition that this copyright notice is kept
%  intact.  V3 by Dave Winkel checks for valid ensembles, skips deviants.


%% Check optional arguments (new: formats(in,out) for multi-filenames, range of seqnos)
if nargin<3
   inrule = {3}; % default 3-digit zero padded file counter;
   % use [] or {0} for single file conversion
end
if nargin<4
   outrule = {1; '.mat'}; % default non-padded out-file counter
end
if nargin<5 | isempty(inrange)
   inrange = [0 inf]; % default, all in-file indices
end
if length(inrange)==1
   inrange = [1 1]*inrange;
end

% Build format arguments for input file names
infmt = ['%s']; inflds = ['transname']; inmult=0;
if iscell(inrule)
   nfi = length(inrule); 
   for i=1:nfi
      x = inrule{i};
      if isnumeric(x) & x>0 % zero-pad for counter
         infmt = [infmt '%0' num2str(x) 'd'];
         inflds = [inflds ',file_no'];
         inmult = inmult+1;
      elseif ischar(x)
         infmt = [infmt '%s'];
         inflds = [inflds ',''' x ''''];
      end
   end
end
if inmult>1
   disp('WARNING: raw2mat_V2, multiple fields for in-file counter');
end
%keyboard

%% Check to see how many input files we have:
nfiles=0; in_files=[];
if inmult<1
   inrange=[1 1]; % single in-file
end
for file_no=inrange(1):inrange(2)
   eval(['inname = sprintf(infmt,' inflds ');']);
   fin=fopen(inname,'r');
   if fin<0
      break;
   end
   fclose(fin);
   nfiles = nfiles+1;
   in_files{nfiles} = inname;
end

if nfiles~=1
   fprintf(1,'Number of in-files found: %d\n',nfiles);
end
if nfiles<1
   return
end

% Build format arguments for output file names
outfmt = ['%s']; outflds = ['matbasename']; outmult=0;
if isempty(outrule)
   outrule={'.mat'};
end
nfo = length(outrule); 
for i=1:nfo
   x = outrule{i};
   if isnumeric(x)
      outfmt = [outfmt '%0' num2str(x) 'd'];
      outflds = [outflds ',matfile'];
      outmult = outmult+1;
   else
      outfmt = [outfmt '%s'];
      outflds = [outflds ',''' x ''''];
   end
end
if outmult>1
   disp('WARNING: raw2mat_V2, multiple fields for out-file counter');
end
%keyboard

% set some parameters now 
MAXENS =900;
MAXBINS = 80; % 128
if outmult<1
   MAXENS = 10000; % 3000; % hope they all fit (FIX later)
end
if nargin<6 | isempty(ens_ind_rng)
    ens_ind_rng = [1 MAXENS];
else
    MAXENS = diff(ens_ind_rng) + 10;
end

% set structures with variable names, their type, and their offset
% within their data block (in bytes)
% I do this so that we have some flexibility should RDI change
% their formats.  It also saves retyping all the names over and over
% during the read and then the write phases of the translation.
name ={'firmwareversion', 'char', 3,
    'firmwarerevison', 'char', 4,
    'sysconfig', 'short', 5,
    'numberbeams','char',9,
    'nbins','char',10,
    'npings','ushort',11,
    'binlen','short',13,
    'blanklen','short',15,
    'watermode','char',17,
    'lowcorr','char',18,
    'codereps','char',19,
    'goodthresh','char',20,
    'errthresh','short',21,
    'tpmin','char',23,
    'tpsec','char',24,
    'tphun','char',25,
    'coords','char',26,
    'headoffset','short',27,
    'headbias','short',29,
    'sensors','char',31,
    'sensorson','char',32,
    'dis1','short',33,
    'pulselen','short',35,
    'fishthresh','char',39,
    'pulselag','short',41,
    'waterband','short',51
};
 offset=name(:,3);type=name(:,2);name=name(:,1);
 fixleader=struct('name',name,'offset',offset,'type',type);
 
 name={
 'ensemble_number', 'ushort', 3,
    'year','char',5,
    'month','char',6,
    'day','char',7,
    'hour','char',8,
    'minute','char',9,
    'second','char',10,
    'hundreths','char',11,
    'ensMSB','char',12,
    'soundspeedRDI','short',15,
    'xducerdepth','short',17,
    'heading','ushort',19,
    'pitch', 'short' ,21,
    'roll', 'short',23,
    'degC', 'short',27,
    'mptmin','char',29,
    'mptsec','char',30,
    'mpthun','char',31,
    'stdhed','char',32,
    'stdpitch','char',33,
    'stdroll','char',34
    };

 offset=name(:,3);type=name(:,2);name=name(:,1);
 varleader=struct('name',name,'offset',offset,'type',type);
 
 name={
 'btnpings','short',3,
  'btdelay','short',5,
  'btcorrthresh','char',7,
  'btminamp','char',8,
  'btpgmin','char',9,
  'btmode','char',10,
  'btmaxerr','short',11,
  'btrange1','ushort',17,
  'btrange2','ushort',19,
  'btrange3','ushort',21,
  'btrange4','ushort',23,
  'btvel1','short',25,
  'btvel2','short',27,
  'btvel3','short',29,
  'btvel4','short',31,
  'btcor1','uchar',33,
  'btcor2','uchar',34,
  'btcor3','uchar',35,
  'btcor4','uchar',36,
  'btevalamp1','char',37,
  'btevalamp2','char',38,
  'btevalamp3','char',39,
  'btevalamp4','char',40,
  'btpg1','char',41,
  'btpg2','char',42,
  'btpg3','char',43,
  'btpg4','char',44,
  'btmaxdep','short',71,
  'rcvstr1','char',73,
  'rcvstr2','char',74,
  'rcvstr3','char',75,
  'rcvstr4','char',76
};
 offset=name(:,3);type=name(:,2);name=name(:,1);
 btleader=struct('name',name,'offset',offset,'type',type);
 
 name={
 'navutcday','uchar',3,
  'navutcmonth','uchar',4,
  'navutcyear','ushort',5,
  'navutctime1','uint',7,
  'navpctimeoff','int',11,
  'navlat1','int',15,
  'navlon1','int',19,
  'navutctime2','uint',23,
  'navlat2','int',27,
  'navlon2','int',31,
  'navavgspd','short',35,
  'navavgtrackT','short',37,
  'navavgtrackM','short',39,
  'navSMG','short',41,
  'navDMG','ushort',43,
  'navflags','ushort',47,
  'navEnsno','uint',51,
  'navadcpyear','ushort',55,
  'navadcpday','uchar',57,
  'navadcpmonth','uchar',58,
  'navadcptime1','uint',59
};  %more later
 offset=name(:,3);type=name(:,2);name=name(:,1);
 navdata=struct('name',name,'offset',offset,'type',type);
 
 % set a bunch of global variables we will need
 trimstr=[];savestr=[];
 for i=1:length(fixleader)
    eval(['global ',fixleader(i).name]);  % make all the fixleader variables global
    % make a string that will let us save these variables whenever we write a 
    % matlab file.
    % Next line added 09-nov-2004, to correct RDI files with 0's in
    % fixedleader of last ensemble:
    trimstr = [trimstr fixleader(i).name '=median(' fixleader(i).name ');'];
    savestr = [savestr ' ' fixleader(i).name];
 end;
 for i=1:length(varleader)
    eval(['global ',varleader(i).name]);
    % make a string that lets us make these variables the right size if
    % we write before the buffer is full...
    trimstr = [trimstr varleader(i).name '=' varleader(i).name '(1:ensembles);'];
    savestr = [savestr ' ' varleader(i).name ];
 end;
%  % Initialize optional Bottom track data and VmDas-(.STA,.LTA) Navg data
%  for i=1:length(btleader)
%     eval(['global ' btleader(i).name]);
%     eval([btleader(i).name ' = [];']);
%     savestr = [savestr ' ' btleader(i).name];
%     trimstr = [trimstr btleader(i).name '=' btleader(i).name ...
%           '(1:min(ensembles,length(' btleader(i).name ' )));'];
%  end
%  for i=1:length(navdata)
%     eval(['global ' navdata(i).name]);
%     eval([navdata(i).name ' = [];']);
%     savestr = [savestr ' ' navdata(i).name];
%     trimstr = [trimstr navdata(i).name '=' navdata(i).name ...
%           '(1:min(ensembles,length(' navdata(i).name ' )));'];
%  end
 
 global v1 v2 v3 v4 e1 e2 e3 e4 cor1 cor2 cor3 cor4 
 global pg1 pg2 pg3 pg4
 % pre-allocate these variables...
 v1=NaN*ones(MAXBINS,MAXENS);v2=v1;v3=v1;v4=v1;
 e1=v1;e2=v1;e3=v1;e4=v1;
 cor1=v1;cor2=v1;cor3=v1;cor4=v1;
 pg1=v1;pg2=v1;pg3=v1;pg4=pg1;
% add them to the savestr...
savestr = [savestr ' v1 v2 v3 v4 e1 e2 e3 e4 cor1 cor2 cor3 cor4 ' ...
    'pg1 pg2 pg3 pg4 '];
% trimstr = [trimstr ' v1=v1(1:nbins,1:ensembles);  v2=v2(1:nbins,1:ensembles);v3= v3(1:nbins,1:ensembles); v4=v4(1:nbins,1:ensembles);' ...
%       'e1=e1(1:nbins,1:ensembles);e2=e2(1:nbins,1:ensembles);e3=e3(1:nbins,1:ensembles);e4=e4(1:nbins,1:ensembles);' ...
%       'cor1=cor1(1:nbins,1:ensembles);cor2=cor2(1:nbins,1:ensembles);cor3=cor3(1:nbins,1:ensembles); cor4=cor4(1:nbins,1:ensembles);'...
%       'pg1=pg1(1:nbins,1:ensembles);pg2=pg2(1:nbins,1:ensembles);pg3=pg3(1:nbins,1:ensembles);pg4=pg4(1:nbins,1:ensembles);' ...
%       'stat1=stat1(1:nbins,1:ensembles);stat2=stat2(1:nbins,1:ensembles);stat3=stat3(1:nbins,1:ensembles); stat4=stat4(1:nbins,1:ensembles);' ...
%       ];
trimstr = [trimstr ' v1=v1(1:nbins,1:ensembles);  v2=v2(1:nbins,1:ensembles);v3= v3(1:nbins,1:ensembles); v4=v4(1:nbins,1:ensembles);' ...
      'e1=e1(1:nbins,1:ensembles);e2=e2(1:nbins,1:ensembles);e3=e3(1:nbins,1:ensembles);e4=e4(1:nbins,1:ensembles);' ...
      'cor1=cor1(1:nbins,1:ensembles);cor2=cor2(1:nbins,1:ensembles);cor3=cor3(1:nbins,1:ensembles); cor4=cor4(1:nbins,1:ensembles);'...
      'pg1=pg1(1:nbins,1:ensembles);pg2=pg2(1:nbins,1:ensembles);pg3=pg3(1:nbins,1:ensembles);pg4=pg4(1:nbins,1:ensembles);' ...
  ];

%% Ready to convert ...
matfile=1;
ensembles=1;
eval(['matname = sprintf(outfmt,' outflds ');']); % first matlab file-name
% now loop through all the files
for file_no=1:nfiles %% For SC files, can really only do one file at a time - dpw 11/04
   % open next in-file
   fprintf(1,'Translating %s\n',in_files{file_no});
   fin=fopen(in_files{file_no},'r','ieee-le'); % open as a little endian
   % Revised to find Ensemble starts=7F7F, where to stop for multi-Ens files.
   % Revised to read large SC files in byte-sized chunks - dpw 11/04
   %%%%%%%%%%%%%%%%%
   iEpos = [];
   ENS_BYTS = 10000; % max expected bytes/RDI_ensem
   RD_NXT = 0; % next read start
   RD_BYTS = 1e6; % read this many bytes for each pass to check for ensembles
   bct = 1;
   while bct>0
       clear A jEpos Enbyt
       ix = fseek(fin, RD_NXT, 'bof');
       if ix<0
           bct = -10; % flag as attempt to move past EOF 
           continue
       end
       [A,bct] = fread(fin, RD_BYTS, 'uchar'); % read another chuck
       %disp('Read'); keyboard
       if bct < 10
           bct = -1;
           continue % not enough more to use (if any)
       end
       % If this chunk went to end-of-file, stop after checking its ensembles:
       if bct < RD_BYTS
           bct = 0; 
       end
       
       jEpos = find(A(1:end-1)==127 & A(2:end)==127); % potential Ens starts
       if ~isempty(jEpos) & jEpos(end)+10 > length(A)
           jEpos(end) = []; % last one won't be complete
       end
       %disp('check'); keyboard
       if isempty(jEpos)
           RD_NXT = RD_NXT + bct - 2; % no valid ensembles in this chunk, move on
           continue
       end
       Enbyt = A(jEpos+3)*256+A(jEpos+2); % bytes in ensemble (w/o checksum)
       ChkAdd=NaN*Enbyt; ChkSum=ChkAdd;
       for iE=1:length(jEpos)
           lb=jEpos(iE)+Enbyt(iE)-1;
           if lb>jEpos(iE) & lb+2<=length(A)
               ChkAdd(iE) = mod(sum(A(jEpos(iE):lb)), 65536); % computed
               ChkSum(iE) = A(lb+2)*256+A(lb+1); % recorded
           end
       end
       iZ=find(ChkAdd-ChkSum==0); % exclude if computed~=recorded
       if isempty(iZ)
           RD_NXT = RD_NXT + bct - 2; % no valid ensembles in this chunk, move on
           continue
       end
       jEpos = jEpos(iZ); Enbyt = Enbyt(iZ);
       % Final check, ALL ensembles should be same length
       EnLEN = median(Enbyt);
       ix = find(Enbyt ~= EnLEN);
       if ~isempty(ix)
           jEpos(ix) = []; Enbyt(ix) = [];
           disp(['   Excluding ' num2str(length(ix)) ' bogus ensemble(s).'])
       end
       % Append ensemble positions, relative to BOF
       if ~isempty(jEpos)
           iEpos = [iEpos; jEpos-1+RD_NXT];
           RD_NXT = iEpos(end) + EnLEN - 2;
           % If last requested ensemble has been found, end search now
           if length(iEpos) > ens_ind_rng(2)
               bct = -2;
           end
       else
           RD_NXT = RD_NXT + bct - 2; % no valid ensembles in this chunk
       end
       clear A
   end
   
   %%%%%%%%%%%%%%%%%
   %keyboard
   Hdrbytes=[]; % Save for diagnosis
   iEb = ens_ind_rng(1); iEe = min(ens_ind_rng(2),length(iEpos));
   for iEns=iEb:iEe % parse specified ensembles - 11/04
      file_pos = iEpos(iEns);
      fseek(fin,file_pos,'bof');
      A=local_hdread(fin);
      ndata=A(1);nbytes=A(2);offsets=A(3:2+ndata);
      % Hdrbytes=[Hdrbytes ndata nbytes offsets];
      %disp(['Hdr ' int2str(iEns) '.']); keyboard
      %if iEns>16 & iEns<22, keyboard; end
      % the raw data is variable length, with variable pieces
      % for storing the data so we need to be flexible in 
      % translating the files.  ndata tells us how many data
      % blocks we have.  offsets tells us where each of these
      % data blocks is in the file stream...
      GotFix = 0; % temporary fix, aeg04 serial VmDas LTAs
      for data_type=1:ndata
         fseek(fin,file_pos+offsets(data_type),'bof');
         data_id=fread(fin,1,'short');
         % disp([iEns data_type offsets(data_type) data_id])
         if ~data_id & GotFix
             disp('Found second fixedleader for this ensemble.')
             keyboard % here (SC), Nav record not really expected
             data_id = 8192; % some Nav records tagged with 00 00
             % in serial binary output; force to be Nav, not fixedleader
         end
         switch data_id
         case 0  % hex bytes = 00 00
            local_fixread(fixleader,fin,file_pos+offsets(data_type),...
                ensembles);
            GotFix = 1;
            %disp(nbins) % two in one ensemble, 2nd=0 after profiles?
         case 128 % 80 00
            local_varread(varleader,fin,file_pos+offsets(data_type),...
               ensembles);
         case 256 % 00 01
            local_velread(fin,file_pos+2+offsets(data_type),...
               ensembles,nbins(ensembles));
         case 512 % 00 02
            local_corread(fin,file_pos+2+offsets(data_type),...
               ensembles,nbins(ensembles));
         case 768 % 00 03
            local_echoread(fin,file_pos+2+offsets(data_type),...
               ensembles,nbins(ensembles));
         case 1024 % 00 04
            local_percread(fin,file_pos+2+offsets(data_type),...
               ensembles,nbins(ensembles));
%          case 1280 % 00 05
%             local_statusread(fin,file_pos+2+offsets(data_type),...
%                ensembles,nbins(ensembles));
%          case 1536 % 00 06
%             local_btread(btleader,fin,file_pos+offsets(data_type),...
%                ensembles);
%          case 8192 % 00 20
%             local_navread(navdata,fin,file_pos+offsets(data_type),...
%                ensembles);
         end; % switch for conditional data reads....
      end; % going through the data types
      % disp(['Ens ct ' int2str(ensembles) ' was read']); keyboard
      % Check to see if ensembles is high enough to write a
      % matlab file yet...
      if ensembles==MAXENS
         eval(trimstr);
         eval(['save ' matname savestr]);
         matfile=matfile+1;
         eval(['matname = sprintf(outfmt,' outflds ');']); % next matlab file-name
         ensembles=0;
      end;
      ensembles=ensembles+1;
      % done with this ensemble, on to next
end; % while going through headers in the file
	%MHA 1/1/99 change: close the file.  If not fid gets updated and it can't do
	%more than 512.
	fclose(fin);
    %keyboard
end; % for file_no=nfiles
% write the last mat file...
%disp('filing...'); keyboard
if ensembles>1 % was 0, DPW 3-2001
   ensembles=ensembles-1;
   eval(trimstr);
   eval(['save ' matname savestr]);
   matfile=matfile+1;
   ensembles=0;
else
    % No ensembles found, just save empty 'ensemble_number'
    ensemble_number=[];
    save(matname, 'ensemble_number');
end;




%%%%%%%%%%%%%%%%%%%%%%%%%% LOCALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ndata]=local_hdread(fin);

offset=ones(1,12)*NaN;
% reads the header information

id = fread(fin,1,'uchar');
if (id~=127)
   ndata=[];
   return;
end;

src= fread(fin,1,'uchar');
nbytes=fread(fin,1,'ushort');
junk=fread(fin,1,'uchar');
ndata=fread(fin,1,'uchar');
for i=1:ndata
   offset(i)=fread(fin,1,'ushort');
end;
ndata = [ndata nbytes offset(1:ndata)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_readfix %%%%%%%%%%%%%%%%%%%%%%
function local_fixread(fixleader,fin,data_pos,ensemble);
% fixleader is set above.  It allows us to define a bunch of variables
% (stored in the name portion of the structure)
for i=1:length(fixleader)
    eval(['global ',fixleader(i).name]);
 end;
 for i=1:length(fixleader)
    fseek(fin,data_pos+fixleader(i).offset-1,'bof');
    
    eval([fixleader(i).name '(ensemble)=fread(fin,1,''' fixleader(i).type ''');']);
    % i.e. firmwareversion=fread(fin,1,'char');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_readvar %%%%%%%%%%%%%%%%%%%%%%
function local_varread(fixleader,fin,data_pos,ensemble);
for i=1:length(fixleader)
    eval(['global ',fixleader(i).name]);
end;
for i=1:length(fixleader)
   fseek(fin,data_pos+fixleader(i).offset-1,'bof');
   eval([fixleader(i).name '(ensemble)=fread(fin,1,''' fixleader(i).type ''');']);
   % i.e. ensemble_number(ensemble)=fread(fin,1,'short');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_velread %%%%%%%%%%%%%%%%%%%%%%
function local_velread(fin,data_pos,ensemble,nbs);
% read in the velocity data... 
if nbs<2, return; end
fseek(fin,data_pos,'bof');
A=reshape(fread(fin,nbs*4,'short'),4,nbs)';
global v1 v2 v3 v4
v1(1:nbs,ensemble)=A(:,1);v2(1:nbs,ensemble)=A(:,2);
v3(1:nbs,ensemble)=A(:,3);v4(1:nbs,ensemble)=A(:,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_velread %%%%%%%%%%%%%%%%%%%%%%
function local_corread(fin,data_pos,ensemble,nbs);
% read in the velocity data... 
if nbs<2, return; end
fseek(fin,data_pos,'bof');
A=reshape(fread(fin,nbs*4,'uchar'),4,nbs)';
global cor1 cor2 cor3 cor4
cor1(1:nbs,ensemble)=A(:,1);cor2(1:nbs,ensemble)=A(:,2);
cor3(1:nbs,ensemble)=A(:,3);cor4(1:nbs,ensemble)=A(:,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_echoread %%%%%%%%%%%%%%%%%%%%%%
function local_echoread(fin,data_pos,ensemble,nbs);
% read in the velocity data... 
if nbs<2, return; end
fseek(fin,data_pos,'bof');
A=reshape(fread(fin,nbs*4,'uchar'),4,nbs)';
global e1 e2 e3 e4
e1(1:nbs,ensemble)=A(:,1);e2(1:nbs,ensemble)=A(:,2);
e3(1:nbs,ensemble)=A(:,3);e4(1:nbs,ensemble)=A(:,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_percread %%%%%%%%%%%%%%%%%%%%%%
function local_percread(fin,data_pos,ensemble,nbs);
% read in the velocity data... 
if nbs<2, return; end
fseek(fin,data_pos,'bof');
A=reshape(fread(fin,nbs*4,'uchar'),4,nbs)';
global pg1 pg2 pg3 pg4
pg1(1:nbs,ensemble)=A(:,1);pg2(1:nbs,ensemble)=A(:,2);
pg3(1:nbs,ensemble)=A(:,3);pg4(1:nbs,ensemble)=A(:,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_statusread %%%%%%%%%%%%%%%%%%%%%%
% function local_statusread(fin,data_pos,ensemble,nbs);
% % read in the velocity data... 
% if nbs<2, return; end
% fseek(fin,data_pos,'bof');
% A=reshape(fread(fin,nbs*4,'uchar'),4,nbs)';
% global stat1 stat2 stat3 stat4
% stat1(1:nbs,ensemble)=A(:,1);stat2(1:nbs,ensemble)=A(:,2);
% stat3(1:nbs,ensemble)=A(:,3);stat4(1:nbs,ensemble)=A(:,4);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_btread %%%%%%%%%%%%%%%%%%%%%%
% function local_btread(fixleader,fin,data_pos,ensemble);
% % same tricks are used here as for local_varread...
% for i=1:length(fixleader)
%     eval(['global ',fixleader(i).name]);
% end;
% for i=1:length(fixleader)
%    fseek(fin,data_pos+fixleader(i).offset-1,'bof');
%    eval([fixleader(i).name '(ensemble)=fread(fin,1,''' fixleader(i).type ''');']);
% end;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% local_navread %%%%%%%%%%%%%%%%%%%%%%
% function local_navread(fixleader,fin,data_pos,ensemble);
% % same tricks are used here as for local_varread...
% for i=1:length(fixleader)
%     eval(['global ',fixleader(i).name]);
% end;
% for i=1:length(fixleader)
%    fseek(fin,data_pos+fixleader(i).offset-1,'bof');
%    eval([fixleader(i).name '(ensemble)=fread(fin,1,''' fixleader(i).type ''');']);
% end;


