function V=ThreeBeam2Earth(jbad,b1,b2,b3,b4,head,pitch,roll,ba,cvxccv,updown)
% Function to convert ADCP beam velocities to Eatth coordinate velocity.
% b1,b2,b3,b4 are four beam velocities
% head is the heading reading, in degrees.
% pitch is the pitch reading, in degrees.
% roll is the roll reading, in degrees.
% ba is the beam angle, in degrees.
% cvxccv is 1 for a convex xdcr, other than 1 for a concave xdcr.  Usually 1.
% updwon is 1 for a down-looking instrument, other than 1 for an up-looker.
%
% Calls matlab routines beam2inst.m and inst2earth.m.
%
i2e=inst2earth(head,pitch,roll,updown);
b2i=beam2inst(ba,cvxccv);
berr = b2i(4,:);
B = [b1 b2 b3 b4]'; 
if jbad == 1;
   B(jbad,:) = -(berr(2)*b2+berr(3)*b3+berr(4)*b4)/berr(1);
end
if jbad == 2;
   B(jbad,:) = -(berr(1)*b1+berr(3)*b3+berr(4)*b4)/berr(2);
end
if jbad == 3;
   B(jbad,:) = -(berr(1)*b1+berr(2)*b2+berr(4)*b4)/berr(3);
end
if jbad == 4;
   B(jbad,:) = -(berr(1)*b1+berr(2)*b2+berr(3)*b3)/berr(4);
end
V = i2e*b2i*B;  % V = [u v w];
return
