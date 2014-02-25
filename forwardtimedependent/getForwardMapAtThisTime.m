%% getForwardMatAtThisTime 
%
% computes the matrix E which gives the forward map at this time that 
% is the concentration at receptor i from source j per unit emission rate 
% is given in E. Therefore we have that 
%
% C(t_k) = E*Q(t_k)
%
% input : 
%           - recept : receptor information and locations
%           - source : source information and locations, note we won't use
%           the emissions here.
%           - wind.dir and wind.vel : velocity and direction of wind.
%           - Vsettling : settling velocity.
%           - stabclass : stability class of atmosphere.

function E = getForwardMapAtThisTime( xdist, ydist, rz, sz, winddir, windvel,...
        Vsettling, Vdeposition, stabclass )
    
    % rotate to go to coordinates of the ermak solution (look at John's
    % paper)
    
    xx = xdist*cos(winddir) + ydist*sin(winddir);
    yy = -xdist*sin(winddir) + ydist*cos(winddir);
    
    % finally compute the matrix E 
    
    E = ermakvec( xx, yy, rz, sz, windvel,Vsettling, Vdeposition, stabclass);
end