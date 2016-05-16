% WHbeam_Process_stg1.m

%% This level of processing makes no corrections for sound speed or clock
%% offsets between the ADCPs and CTD (or GPS).  Time is corrected, however,
%% for jumps in the ADCP clock noted in the matlab-datafile index
%% (i.e., ADDN jumping ahead 1 hour in home02 and bs03).
%% Pitch and roll attitude of the ADCPs is factored in to the velocity

% initialize output structure
clear ADP
ADP.u = NaN * ADRaw.v1; ADP.v = ADP.u; ADP.w = ADP.u; ADP.werr = ADP.u;
ADP.ec1 = 0*ones(size(ADP.u)); ADP.ec2 = ADP.ec1;
ADP.ec3 = ADP.ec1; ADP.ec4 = ADP.ec1;
ADP.yday = ADRaw.yday;
ADP.ens_no = ADRaw.ens_no;
ADP.soundvel = ADRaw.soundvel;
ADP.svel_calc = ADP.soundvel;
ADP.z_adcp = ADRaw.z_adcp;
ADP.depth_xducer = ADRaw.depth_xducer;
vvs = {'svel_calc','ens_up','yd_up','temp'};
for i=1:length(vvs)
    if isfield(ADRaw, vvs{i}) && ~isempty(ADRaw.(vvs{i}))
        ADP.(vvs{i}) = ADRaw.(vvs{i});
    end
end
if size(ADP.z_adcp,1)==1
    ADP.z_adcp = ADP.z_adcp'; % column vector
end

% Screen data for 'gross' problems
ix = find(ADRaw.cor1_bm(:) < WC_val);
ADRaw.v1(ix) = NaN;
ix = find(ADRaw.cor2_bm(:) < WC_val);
ADRaw.v2(ix) = NaN;
ix = find(ADRaw.cor3_bm(:) < WC_val);
ADRaw.v3(ix) = NaN;
ix = find(ADRaw.cor4_bm(:) < WC_val);
ADRaw.v4(ix) = NaN;

% RDI transformation matrix from beams 1-4 to
%   u(1-2), v(4-3), w(avg xz,yz), err vel 
Bm2InTx = beam2inst(theta_o, Cnvx);
% Make range vector from the depth vector (before pitch/roll adjustments)
c_tho = cos(theta_o*pi/180); s_tho = sin(theta_o*pi/180);
rBM_o = ADP.z_adcp/c_tho;

% % Turn the recorded pitch into real pitch... insignif. for small rolls.
% ADRaw.pitch=180/pi*atan(tan(pi/180*ADRaw.pitch).*cos(pi/180*ADRaw.roll));
% % Slighty filter pitch and roll:
% ADP.pitch = [ ADRaw.pitch(1), ...
%         (0.2*ADRaw.pitch(1:end-2)) + (0.6*ADRaw.pitch(2:end-1)) + (0.2*ADRaw.pitch(3:end)), ...
%         ADRaw.pitch(end) ];
% ADP.roll = [ ADRaw.roll(1), ...
%         (0.2*ADRaw.roll(1:end-2)) + (0.6*ADRaw.roll(2:end-1)) + (0.2*ADRaw.roll(3:end)), ...
%         ADRaw.roll(end) ];
%% Just copy pitch and roll (and compass heading) instead:
ADP.pitch = ADRaw.pitch;
ADP.roll = ADRaw.roll;
ADP.heading = ADRaw.heading;
% Signs for correcting range(depth)
if UpDown == 1
    sg1 = 1; sg3 = 1;
else
    sg1 = 1; sg3 = -1;
end

%% Loop thru pings (ensembles), do bin mapping prior transforming beams-to-earth
warning off MATLAB:interp1:NaNinY
for ic=1:length(ADP.yday)
    % The first step is to compute depth vectors for each beam.  These differ from each other 
    % since the beams are tilted.  DO NOT adjust for sound speed (this will happen later).
    % Roll affects beams 1 and 2,  pitch affects beams 3 and 4.
    rBM = rBM_o; % adjust depths for pitched/rolled beam angles
    z1_rel = cos( pi/180*(theta_o + sg1*ADP.roll(ic)) )*rBM;
    z2_rel = cos( pi/180*(theta_o + -sg1*ADP.roll(ic)) )*rBM; 
    z3_rel = cos( pi/180*(theta_o + sg3*ADP.pitch(ic)) )*rBM;
    z4_rel = cos( pi/180*(theta_o + -sg3*ADP.pitch(ic)) )*rBM;
    
    % map beam velocities onto standard depths (extend first bin nearer, if needed)
    u1 = interp1( [0; z1_rel], [ADRaw.v1(1,ic); ADRaw.v1(:,ic)], ADP.z_adcp);
    u2 = interp1( [0; z2_rel], [ADRaw.v2(1,ic); ADRaw.v2(:,ic)], ADP.z_adcp);
    u3 = interp1( [0; z3_rel], [ADRaw.v3(1,ic); ADRaw.v3(:,ic)], ADP.z_adcp);
    u4 = interp1( [0; z4_rel], [ADRaw.v4(1,ic); ADRaw.v4(:,ic)], ADP.z_adcp);
    % map echo intensities, for later screening (extend first bin nearer, if needed)
    ADP.ec1(:,ic) = round( interp1( [0; z1_rel], [ADRaw.ec1_bm(1,ic); ADRaw.ec1_bm(:,ic)], ...
        ADP.z_adcp, 'linear', 0) );
    ADP.ec2(:,ic) = round( interp1( [0; z2_rel], [ADRaw.ec2_bm(1,ic); ADRaw.ec2_bm(:,ic)], ...
        ADP.z_adcp, 'linear', 0) );
    ADP.ec3(:,ic) = round( interp1( [0; z3_rel], [ADRaw.ec3_bm(1,ic); ADRaw.ec3_bm(:,ic)], ...
        ADP.z_adcp, 'linear', 0) );
    ADP.ec4(:,ic) = round( interp1( [0; z4_rel], [ADRaw.ec4_bm(1,ic); ADRaw.ec4_bm(:,ic)], ...
        ADP.z_adcp, 'linear', 0) );
    
    % Now, transform beam velocities to earth (geomagnetic) coordinates.
    In2Geo = inst2earth(ADP.heading(ic), ADP.pitch(ic), ADP.roll(ic), UpDown);
    VelsGeo = ( In2Geo * Bm2InTx * [u1 u2 u3 u4]' )';
    % Gather components into output structure
    ADP.u(:,ic) = VelsGeo(:,1);
    ADP.v(:,ic) = VelsGeo(:,2);
    ADP.w(:,ic) = VelsGeo(:,3);
    ADP.werr(:,ic) = VelsGeo(:,4);
    
end % of transforming pings from beam to earth coordinates

%clear ADRaw
% %% Trim to original yearday range
% vn = fieldnames(ADP);
% id = find( ADP.yday>=yd_b & ADP.yday<=yd_e );
% for iv=1:length(vn)
%     if ~strcmp(vn{iv}, 'z_adcp')
%         ADP.(vn{iv}) = ADP.(vn{iv})(:,id);
%     end
% end

return

%% Code to compute without Cor threshold, bin mapping:
clear vRDI
for ic=1:length(ADRaw.yday)
    vRDI = Beam2Earth(ADRaw.v1(:,ic),ADRaw.v2(:,ic),ADRaw.v3(:,ic),ADRaw.v4(:,ic),ADRaw.heading(ic),ADRaw.pitch(ic),ADRaw.roll(ic), 20,1,1);
    vRD.u(:,ic) = vRDI(1,:)';
    vRD.v(:,ic) = vRDI(2,:)';
    vRD.w(:,ic) = vRDI(3,:)';
    vRD.werr(:,ic) = vRDI(4,:)';
    clear vRDI
end

%%
if 0 %% copy and paste to compute smoothed vels:
    npts = 180;
    for i=1:length(ADP.z_adcp)
        SM.uRD(i,:) = medfilt1(vRD.u(i,:), npts);
        SM.vRD(i,:) = medfilt1(vRD.v(i,:), npts);
        SM.wRD(i,:) = medfilt1(vRD.w(i,:), npts);
        SM.uDW(i,:) = medfilt1(ADP.u(i,:), npts);
        SM.vDW(i,:) = medfilt1(ADP.v(i,:), npts);
        SM.wDW(i,:) = medfilt1(ADP.w(i,:), npts);
    end
end