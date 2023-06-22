function data = convert_units_SBE21(rawData,config)
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
T = ITS90*S + O;
%
% Convert conductivity: Use Coeffs1 bce flag UseG_J==1.
G = sscanf(config.ConductivitySensor.Coefficients1.G,'%e');
H = sscanf(config.ConductivitySensor.Coefficients1.H,'%e');
I = sscanf(config.ConductivitySensor.Coefficients1.I,'%e');
J = sscanf(config.ConductivitySensor.Coefficients1.J,'%e');
CPcor = sscanf(config.ConductivitySensor.Coefficients1.CPcor,'%e');
CTcor = sscanf(config.ConductivitySensor.Coefficients1.CTcor,'%e');
%
C = (G + H*rawData.cond.^2 + I*rawData.cond.^3 + J*rawData.cond.^4)./(10*(1 + CTcor*T + CPcor*rawData.pres));
% Auxilliary Sensors:
%
data = struct('T',T,'C',C);%,'S',S,'AFL',AFL,'FTU',FTU,'Flo',Flo);
end