%test code
cd /Volumes/MHA-STICK/alarms/

%initialize
handles.bthits=0;
handles.btrange=nan;
handles.decreasing=0;

handles.max_warn_freq_sec= 60; %max frequency in secs of warnings.
handles.hits_to_warn=3; %required successive bottom detections to warn

handles.last=now - handles.max_warn_freq_sec/3600/24; %set the time of the last bottom warning in the past

%% Generate a test loop of bottom ranges simulating a bottom approach
vals=[nan 108 104 100 96 nan 100 96 92 60 30 nan 34 36 90 100 ]*100;
vals=[nan(1,3) 140:-4:30 35:4:100 nan];
vals(20)=nan;

for c=1:length(vals)
    if handles.enable_bottom_warnings
        
        Vel.btrange=vals(c)
        
        %%
        %Vel.btrange=nan;
        
        %%
        disp(['range is ' num2str(vals(c))])
        thisrange=nanmean(Vel.btrange)/100;
        if ~isnan(thisrange)  %if last ping we saw the bottom
            if thisrange < handles.btrange
                handles.decreasing=1
            else
                handles.decreasing=0;
            end
            
            handles.btrange= thisrange; %set the new bottom range - average of all four beams, in m
            handles.bthits=handles.bthits+1
            disp 'got one'
        else
            handles.decreasing=0;
            handles.bthits=0
            disp 'skipping.'
        end
        
        if handles.bthits >=handles.hits_to_warn & handles.decreasing & now > handles.last + handles.max_warn_freq_sec/3600/24
            %get people's attention...
            %%
            mrt=imread('MrT.jpg');
            Y4=wavread(['mistake_hg.wav']);
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
    end %end if bottom warnings enabled
    disp([num2str((now - handles.last)*24*3600) ' s since last message.'])
    pause(3)
    
end
%& thisrange < handles.btrange