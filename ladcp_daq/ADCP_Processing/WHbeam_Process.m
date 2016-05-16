% WHbeam_Process.m

% initialize output structure
clear ADP
ADP.depth = ZGrid; % surface-relative depth grid
Trow = length(ADP.depth); Tcol = length(ADRaw.dtnum);
ADP.u = NaN * ones(Trow,Tcol); ADP.v = ADP.u;
ADP.w = ADP.u; ADP.werr = ADP.u;
ADP.ec1 = 0*ones(Trow,Tcol); ADP.ec2 = ADP.ec1;
ADP.ec3 = ADP.ec1; ADP.ec4 = ADP.ec1;
ADP.dtnum = ADRaw.dtnum;
ADP.ens_no = ADRaw.ens_no;
ADP.soundvel = ADRaw.soundvel;
ADP.svel_calc = ADP.soundvel;

ADP.z_adcp = -ADRaw.z_adcp;  %sign changes depending on up or down-looking. (negative down, positive up...need to fix).
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
% soundspeed correction, actual/nominal
if isempty(ADP.svel_calc)
    ADP.svel_calc = ADP.soundvel;
end
SSadj = ADP.svel_calc ./ ADP.soundvel; 
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

%% Loop thru pings (ensembles), do bin mapping before transforming beams-to-earth
warning off MATLAB:interp1:NaNinY
for ic=1:length(ADP.dtnum)
    % The first step is to compute depth vectors for each beam.  These differ from 
    % each other since the beams are tilted; Then adjust for sound speed.
    % Roll affects beams 1 and 2,  pitch affects beams 3 and 4.
    rBM = rBM_o; % adjust depths for pitched/rolled beam angles
    z1_rel = cos( pi/180*(theta_o + sg1*ADP.roll(ic)) )*rBM;
    z2_rel = cos( pi/180*(theta_o + -sg1*ADP.roll(ic)) )*rBM; 
    z3_rel = cos( pi/180*(theta_o + sg3*ADP.pitch(ic)) )*rBM;
    z4_rel = cos( pi/180*(theta_o + -sg3*ADP.pitch(ic)) )*rBM;
    % compute absolute depths considering sound speed and ADCP depth
    z1_abs = z1_rel; z1_abs(1) = z1_rel(1) * SSadj(ic) + ADP.depth_xducer(ic);
    z1_abs(2:end) = z1_abs(1) + cumsum([diff(z1_rel)*SSadj(ic)]);
    z2_abs = z2_rel; z2_abs(1) = z2_rel(1) * SSadj(ic) + ADP.depth_xducer(ic);
    z2_abs(2:end) = z2_abs(1) + cumsum([diff(z2_rel)*SSadj(ic)]);
    z3_abs = z3_rel; z3_abs(1) = z3_rel(1) * SSadj(ic) + ADP.depth_xducer(ic);
    z3_abs(2:end) = z3_abs(1) + cumsum([diff(z3_rel)*SSadj(ic)]);
    z4_abs = z4_rel; z4_abs(1) = z4_rel(1) * SSadj(ic) + ADP.depth_xducer(ic);
    z4_abs(2:end) = z4_abs(1) + cumsum([diff(z4_rel)*SSadj(ic)]);
    % map beam velocities onto standard depths (sfc-refn'd),
    %   adjust for actual soundspeed
    u1 = interp1( z1_abs, ADRaw.v1(:,ic), ADP.depth) * SSadj(ic);
    u2 = interp1( z2_abs, ADRaw.v2(:,ic), ADP.depth) * SSadj(ic);
    u3 = interp1( z3_abs, ADRaw.v3(:,ic), ADP.depth) * SSadj(ic);
    u4 = interp1( z4_abs, ADRaw.v4(:,ic), ADP.depth) * SSadj(ic);
    % map echo intensities, for later screening
    ADP.ec1(:,ic) = round( interp1( z1_abs, ADRaw.ec1_bm(:,ic), ...
        ADP.depth, 'linear', 0) );
    ADP.ec2(:,ic) = round( interp1( z2_abs, ADRaw.ec2_bm(:,ic), ...
        ADP.depth, 'linear', 0) );
    ADP.ec3(:,ic) = round( interp1( z3_abs, ADRaw.ec3_bm(:,ic), ...
        ADP.depth, 'linear', 0) );
    ADP.ec4(:,ic) = round( interp1( z4_abs, ADRaw.ec4_bm(:,ic), ...
        ADP.depth, 'linear', 0) );
    
    % Now, transform beam velocities to earth (geomagnetic) coordinates.
    In2Geo = inst2earth(ADP.heading(ic)+MagDec, ADP.pitch(ic), ADP.roll(ic), UpDown);
    VelsGeo = ( In2Geo * Bm2InTx * [u1 u2 u3 u4]' )';
    % Gather components into output structure
    ADP.u(:,ic) = VelsGeo(:,1);
    ADP.v(:,ic) = VelsGeo(:,2);
    ADP.w(:,ic) = VelsGeo(:,3);
    ADP.werr(:,ic) = VelsGeo(:,4);
    
end % of transforming pings from beam to earth coordinates

return
