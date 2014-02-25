%% getDepositionAtReceptors

function d = getDepositionAtDustfallJar( dustfalljar, A, dt, M, Vd )

d = (A * dt * M * Vd)*[dustfalljar(:).accum_concentration]';

end