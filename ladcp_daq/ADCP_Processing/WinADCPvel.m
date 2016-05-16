   clear
   clf 
   bt = datenum(2008,6,13,8,50,0);
   et = datenum(2008,6,13,9,50,0);
%  bt = datenum(2008,6,13,00,00,0);
%  et = datenum(2008,6,14,06,00,0);
%  bt = datenum(2008,6,13,17,10,0);
%  et = datenum(2008,6,13,17,20,0);
%  bt = datenum(2008,6,11,20,0,0);
%  et = datenum(2008,6,11,20,30,0);
%  bt = datenum(2008,6,12,8,10,0);
%  et = datenum(2008,6,12,8,40,0);
%  bt = datenum(2008,6,13,12,50,0);
%  et = datenum(2008,6,13,13,40,0);
   lpdt = 120; plot_300=1; plot_1200=0;
 A=load('../winadcp/AD3dn_mat/OR3DN');A.dz=1;A.nz=33;A.dt=0.5;A.z0=15;
 B=load('../winadcp/AD12up_mat/OR12UP');B.dz=-0.5;B.nz=40;B.dt=0.33;B.z0=29.75;
 C=load('../winadcp/AD12dn_mat/OR12DN');C.dz=0.5;C.nz=34;C.dt=0.33;C.z0=31;

   A.z = A.z0+(1:A.nz)'*A.dz;
   A.pitch = A.AnP100thDeg/100; A.roll = A.AnR100thDeg/100;  
   A.heading = A.AnT100thDeg/100; A.u=A.SerEmmpersec'/1000; 
   A.v=A.SerNmmpersec'/1000; A.w=A.SerVmmpersec'/1000; 
   A.ver=A.SerErmmpersec'/1000; A.Pr = A.AnDepthmm/1000;
   A.Jday_lct = datenum(2008,6,A.SerDay,A.SerHour,A.SerMin,A.SerSec)-7/24;
   A.jj = find(A.Jday_lct >=bt & A.Jday_lct <= et);
   for i = 1:A.nz;
       A.lpu = lpass(A.u(i,A.jj),A.dt,lpdt,4,'b'); A.lpU(i,:) = A.lpu(:)';  
       A.lpv = lpass(A.v(i,A.jj),A.dt,lpdt,4,'b'); A.lpV(i,:) = A.lpv(:)'; 
       A.lpw = lpass(A.w(i,A.jj),A.dt,lpdt,4,'b'); A.lpW(i,:) = A.lpw(:)'; 
   end
%  A.ii = find(A.Jday_lct(A.jj)>=(bt+1/60/24) & A.Jday_lct(A.jj)<=(bt+2/60/24));
%  A.lpU = A.lpU-mymean(A.lpU(:,A.ii)')'*ones(1,size(A.lpU,2));
%  A.lpV = A.lpV-mymean(A.lpV(:,A.ii)')'*ones(1,size(A.lpV,2));
%  A.lpW = A.lpW-mymean(A.lpW(:,A.ii)')'*ones(1,size(A.lpW,2));

   B.z = B.z0+(1:B.nz)'*B.dz;
   B.pitch = B.AnP100thDeg/100; B.roll = B.AnR100thDeg/100;  
   B.heading = B.AnT100thDeg/100; B.u=B.SerEmmpersec'/1000; 
   B.v=B.SerNmmpersec'/1000; B.w=B.SerVmmpersec'/1000; 
   B.ver=B.SerErmmpersec'/1000; B.Pr = B.AnDepthmm/1000;
   B.Jday_lct = datenum(2008,6,B.SerDay,B.SerHour,B.SerMin,B.SerSec)-7/24;
   B.jj = find(B.Jday_lct >=bt & B.Jday_lct <= et);
   for i = 1:B.nz;
       B.lpu = lpass(B.u(i,B.jj),B.dt,lpdt,4,'b'); B.lpU(i,:) = B.lpu(:)';  
       B.lpv = lpass(B.v(i,B.jj),B.dt,lpdt,4,'b'); B.lpV(i,:) = B.lpv(:)'; 
       B.lpw = lpass(B.w(i,B.jj),B.dt,lpdt,4,'b'); B.lpW(i,:) = B.lpw(:)'; 
   end
%  B.ii = find(B.Jday_lct(B.jj)>=(bt+1/60/24) & B.Jday_lct(B.jj)<=(bt+2/60/24));
%  B.lpU = B.lpU-mymean(B.lpU(:,B.ii)')'*ones(1,size(B.lpU,2));
%  B.lpV = B.lpV-mymean(B.lpV(:,B.ii)')'*ones(1,size(B.lpV,2));
%  B.lpW = B.lpW-mymean(B.lpW(:,B.ii)')'*ones(1,size(B.lpW,2));

   C.z = C.z0+(1:C.nz)'*C.dz;
   C.pitch = C.AnP100thDeg/100; C.roll = C.AnR100thDeg/100;  
   C.heading = C.AnT100thDeg/100; C.u=C.SerEmmpersec'/1000; 
   C.v=C.SerNmmpersec'/1000; C.w=C.SerVmmpersec'/1000; 
   C.ver=C.SerErmmpersec'/1000; C.Pr = C.AnDepthmm/1000;
   C.Jday_lct = datenum(2008,6,C.SerDay,C.SerHour,C.SerMin,C.SerSec)-7/24;
   C.jj = find(C.Jday_lct >=bt & C.Jday_lct <= et);
   for i = 1:C.nz;
       C.lpu = lpass(C.u(i,C.jj),C.dt,lpdt,4,'b'); C.lpU(i,:) = C.lpu(:)';  
       C.lpv = lpass(C.v(i,C.jj),C.dt,lpdt,4,'b'); C.lpV(i,:) = C.lpv(:)'; 
       C.lpw = lpass(C.w(i,C.jj),C.dt,lpdt,4,'b'); C.lpW(i,:) = C.lpw(:)'; 
   end
%  C.ii = find(C.Jday_lct(C.jj)>=(bt+1/60/24) & C.Jday_lct(C.jj)<=(bt+2/60/24));
%  C.lpU = C.lpU-mymean(C.lpU(:,C.ii)')'*ones(1,size(C.lpU,2));
%  C.lpV = C.lpV-mymean(C.lpV(:,C.ii)')'*ones(1,size(C.lpV,2));
%  C.lpW = C.lpW-mymean(C.lpW(:,C.ii)')'*ones(1,size(C.lpW,2));
   npt = min([size(B.lpU,2) size(C.lpU,2)]);
   if length(B.jj) >= length(C.jj); Jday_lct= C.Jday_lct(C.jj); end
   if length(B.jj) <= length(C.jj); Jday_lct= B.Jday_lct(B.jj); end

   clf
   x0 = 0.1; y0 = 0.7; dx = 0.7; dy = 0.23; ddy = 0.01;
   ax = axes('position',[x0 y0 dx dy],'ydir','reverse','box','on');
   hold on
   uc = -0.6:0.04:0.2; vc = -0.6:0.04:0.2; wc = -0.2:0.02:0.2;
   lpU = [flipud(B.lpU(:,1:npt));C.lpU(:,1:npt)];
   z = [flipud(B.z(:)); C.z(:)];
   if plot_300; [c,h]=contourf(A.Jday_lct(A.jj),A.z,A.lpU,uc);  end
%  [c,h]=contourf(B.Jday_lct(B.jj),B.z,B.lpU,uc); 
%  [c,h]=contourf(C.Jday_lct(C.jj),C.z,C.lpU,uc); 
   if plot_1200; [c,h]=contourf(Jday_lct,z,lpU,uc);  end
   caxis([uc(1) uc(end)])
   myaxis;
   datetick('x','keepticks','keeplimits');
   hc = colorbar('vertical');
   set(hc,'position',[x0+dx+0.01 y0 0.02 dy]);
   text(1.12,0.3,'U (m s^{-1})','unit','normalized',...
             'rotation',90);
   set(ax,'position',[x0 y0 dx dy]);
   ylabel('Depth (m)');
   title(datestr(bt,2));

   y0 = y0-dy-ddy;
   ax = axes('position',[x0 y0 dx dy],'ydir','reverse','box','on');
   hold on
   lpV = [flipud(B.lpV(:,1:npt));C.lpV(:,1:npt)];
   if plot_300; [c,h]=contourf(A.Jday_lct(A.jj),A.z,A.lpV,vc); end
%  [c,h]=contourf(B.Jday_lct(B.jj),B.z,B.lpV,vc); 
%  [c,h]=contourf(C.Jday_lct(C.jj),C.z,C.lpV,vc); 
   if plot_1200; [c,h]=contourf(Jday_lct,z,lpV,vc); end
   caxis([vc(1) vc(end)])
   myaxis;
   datetick('x','keepticks','keeplimits');
   hc = colorbar('vertical');
   set(hc,'position',[x0+dx+0.01 y0 0.02 dy]);
   text(1.12,0.3,'V (m s^{-1})','unit','normalized',...
             'rotation',90);
   set(ax,'position',[x0 y0 dx dy]);
   ylabel('Depth (m)');

   y0 = y0-dy-ddy;
   ax = axes('position',[x0 y0 dx dy],'ydir','reverse','box','on');
   hold on
   lpW = [flipud(B.lpW(:,1:npt));C.lpW(:,1:npt)];
   if plot_300; [c,h]=contourf(A.Jday_lct(A.jj),A.z,A.lpW,wc);  end
%  [c,h]=contourf(B.Jday_lct(B.jj),B.z,B.lpW,wc); 
%  [c,h]=contourf(C.Jday_lct(C.jj),C.z,C.lpW,wc); 
   if plot_1200; [c,h]=contourf(Jday_lct,z,lpW,wc); end
   caxis([wc(1) wc(end)])
   myaxis;
   datetick('x','keepticks','keeplimits');
   hc = colorbar('vertical');
   set(hc,'position',[x0+dx+0.01 y0 0.02 dy]);
   text(1.12,0.3,'W (m s^{-1})','unit','normalized',...
             'rotation',90);
   set(ax,'position',[x0 y0 dx dy]);
   ylabel('Depth (m)');
