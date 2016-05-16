function V=Beam2Earth(b1,b2,b3,b4,head,pitch,roll,ba,cvxccv,updown)
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
V = i2e*b2i*[b1 b2 b3 b4]'; %V = [u v w V_err];
return
