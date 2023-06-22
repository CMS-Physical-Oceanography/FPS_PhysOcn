function data = convert_units_SBE25(rawData,config)



% Convert temperature to ITS-90:
% 1) extract calibration coefficients
G = sscanf(config.TemperatureSensor.G,'%e');
H = sscanf(config.TemperatureSensor.H,'%e');
I = sscanf(config.TemperatureSensor.I,'%e');
J = sscanf(config.TemperatureSensor.J,'%e');
F0= sscanf(config.TemperatureSensor.F0,'%e');
S = sscanf(config.TemperatureSensor.Slope,'%f');
O = sscanf(config.TemperatureSensor.Offset,'%f');
% 2) evaluate calibration equation
ITS90 = 1./(G+H*log(F0./rawData.temp)+I*log(F0./rawData.temp).^2+J*log(F0./rawData.temp).^3) - 273.15;
% 3) apply slope/offset correction
ITS90 = ITS90*S + O;
%
% Convert pressure to (lbs/in^2 absolute; PSIA)
% 1) extract calibration coefficients
PA0     = sscanf(config.PressureSensor.PA0,'%e');
PA1     = sscanf(config.PressureSensor.PA1,'%e');
PA2     = sscanf(config.PressureSensor.PA2,'%e');
PTEMPA0 = sscanf(config.PressureSensor.PTEMPA0,'%e');
PTEMPA1 = sscanf(config.PressureSensor.PTEMPA1,'%e');
PTEMPA2 = sscanf(config.PressureSensor.PTEMPA2,'%e');
PTCA0   = sscanf(config.PressureSensor.PTCA0,'%e');
PTCA1   = sscanf(config.PressureSensor.PTCA1,'%e');
PTCA2   = sscanf(config.PressureSensor.PTCA2,'%e');
PTCB0   = sscanf(config.PressureSensor.PTCB0,'%e');
PTCB1   = sscanf(config.PressureSensor.PTCB1,'%e');
PTCB2   = sscanf(config.PressureSensor.PTCB2,'%e');
O       = sscanf(config.PressureSensor.Offset,'%e');
% 2) Evaluate calibration equations;
y    = rawData.ptmp;
t    = PTEMPA0 + PTEMPA1*y + PTEMPA2*y.^2;
x    = rawData.pres - (PTCA0 + PTCA1*t + PTCA2*t.^2);
n    = PTCB0.*x./(PTCB0 + PTCB1*t + PTCB2*t.^2);
PSIA = PA0 + PA1*n + PA2*n.^2;
% 3) Convert to Pascals, then decibars
P    = PSIA*(6894.76./1e4);
%


% Convert conductivity to psu



% Convert dissolved oxygen
% 0) Oxygen Solubility based on Benson & Krause's refit of Garcia & Gordon 1992 
Ts = ln( (298.15-ITS90)./(273.15+ITS90) );
A0 = 2.00907; B0 = -0.00624523;
A1 = 3.22014; B1 = -0.00737614;
A2 = 4.0501;  B2 = -0.0103410;
A3 = 4.94457; B3 = -0.00817083
A4 = -0.256847;
A5 = 3.88767;
C0 = -0.00000048868;
Oxsol = exp( A0 + A1*Ts + A2*Ts.^2 + A3*Ts.^3 + A4*Ts.^4 + A5*Ts.^5 +...
             S.*( B0 + B1*Ts + B2*Ts.^2 + B3*Ts.^3 ) + C0*S.^2 );
% 1) allocate coefficients from con file
Soc = sscanf(config.OxygenSensor.Coefficients1.Soc,'%e');
A   = sscanf(config.OxygenSensor.Coefficients1.A,'%e');
B   = sscanf(config.OxygenSensor.Coefficients1.B,'%e');
C   = sscanf(config.OxygenSensor.Coefficients1.C,'%e');
D0  = sscanf(config.OxygenSensor.Coefficients1.D0,'%e');
D1  = sscanf(config.OxygenSensor.Coefficients1.D1,'%e');
D2  = sscanf(config.OxygenSensor.Coefficients1.D2,'%e');
E   = sscanf(config.OxygenSensor.Coefficients1.E,'%e');
H1  = sscanf(config.OxygenSensor.Coefficients1.H1,'%e');
H2  = sscanf(config.OxygenSensor.Coefficients1.H2,'%e');
H3  = sscanf(config.OxygenSensor.Coefficients1.H3,'%e');
Tau20  = sscanf(config.OxygenSensor.Coefficients1.Tau20,'%e');
offset = sscanf(config.OxygenSensor.Coefficients1.offset,'%e');
% 2) evaluate
DO2 = Soc*(rawData.aux0+offset).*(1.0 + A*T + B*T.^2 + C*T.^3).*Oxsol.*exp(E*P./(T+273.15));




data = struct('T',T,'C',C,'S',S,'P',P,'dO',dO,'PAR',PAR,'AFL',AFL,'FTU',FTU,'Flo',Flo);
