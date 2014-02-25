%% test timedepsource against previous code
%
% 

% call forward solver
%
% make sure source.Q is set to be constant in forwardtimedep script.
forwardtimedep;

% now use old code to compute depositions and compare solutions 

% this is silly but the code is setup to work with row vectors
source.x = source.x';
source.y = source.y'; 
source.z = source.z';
source.Q = ones(1, source.n);
recept.x = recept.x';
recept.y = recept.y';
recept.z = recept.z';

dep = zeros(source.n, recept.n);
for k = [1 : tskip : nwind ]
   
    Vw = wind.vel(k);
    theta = wind.dir(k);
    
    dep = getDepAtRecept(dep, recept, source, theta, Vw, Vspbo, Vdpbo, ...
        stabclass, A, Mpb, dt);
    
    
end

% sum over deposition from all source
dep = sum(dep, 1);

figure
bar(dep*1e6);
xlabel('receptor');
ylabel('deposition in mg');