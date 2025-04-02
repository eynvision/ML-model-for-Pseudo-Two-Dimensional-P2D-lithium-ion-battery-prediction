function potential=equi_LiMnO(x)
global PM

potential=interp1(PM.LiMnO_data(1,:),PM.LiMnO_data(2,:),real(x));
end