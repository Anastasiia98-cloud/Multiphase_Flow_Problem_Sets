function [c] = updateSoundSpeed(p,rho,T,Yv,Yl)

% Definition of parameters from the SG parameter vector
[Cv,Cl,gammav,gammal,~,~,pinfv,pinfl,~,~] = stiffenedGasParameters;

% Calculation of gas and liquid densities and volume fractions
rhov = (p+pinfv)/((gammav-1)*Cv*T);
rhol = (p+pinfl)/((gammal-1)*Cl*T);
alphav = Yv*rho/rhov;
alphal = Yl*rho/rhol;

cv_squared = gammav*(p+pinfv)/rhov; % Squared speed of sound of the single gas phase
cl_squared = gammal*(p+pinfl)/rhol; % Squared speed of sound of the single liquid phase
c_Wood = sqrt(1/(rho*(alphal/(rhol*cl_squared)+alphav/(rhov*cv_squared)))); % Wood speed of sound

% Exact mixture speed of sound of 4-equation model
c = 1/sqrt(1/c_Wood^2+rho*T*alphav*alphal*rhov*rhol*gammav*gammal*Cv*Cl/(alphav*rhov*gammav*Cv+alphal*rhol*gammal*Cl)*((gammav-1)/(rhov*cv_squared)+(gammal-1)/(rhol*cl_squared))^2);

end