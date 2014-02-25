%% getDepositionAtReceptors

function d = getDepositionAtReceptors( dustfalljar, A, dt, M, Vd )

d = (A * dt * M * Vd)*[dustfalljar(:).accum_concentration]';

end