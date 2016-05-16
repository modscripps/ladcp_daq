function y=Earth2Beam(ve,vn,vu,verr,head,pitch,roll,ba,cvxccv,updown)
% Function to convert ADCP earth coordinate velocities to beam velocity.
% ve is East Velocity, vn is North Velocity, vu is Up Velocity, verr is Error Velocity
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
e2i=inv(i2e);
b2i=beam2inst(ba,cvxccv);
i2b=inv(b2i);
y=i2b*e2i*[ve vn vu verr]';
return
