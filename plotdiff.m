

model.tmod = 'V2F';
model.fmod = {'DWX'};
model.cases= {'H2O','CO2','H2OP','H2O','CO2','H2OP','H2O','CO2','H2OP'};
model.rmod = [zeros(3,1); ones(6,1)];
model.kmod = [zeros(6,1); 2*ones(3,1)];
model.dens = ones(11,1);


plotcases('t10r',{'V2F'},{'DWX'},[0 1],[0])