function potential=equi_LixC(x)
potential=8.00229+5.0647*x-12.578*x.^0.5-8.6322e-4*x.^-1+2.1765e-5*x.^1.5-0.46016*exp(15*(0.06-x))-0.55364*exp(-2.4326*(x-0.92));
end