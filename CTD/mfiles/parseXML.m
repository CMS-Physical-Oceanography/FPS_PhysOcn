function config = parseXML(conFile)
% PARSEXML Convert XML file to a MATLAB structure.

fid = fopen(conFile);

levl = 0;
note = 0;
coeff= 0;
config = struct([]);
while ~feof(fid)
    tline = fgetl(fid);
    if levl==0
        sensorStr = regexp(tline, '(?<=<Sensor\s).*(?=>)', 'match');
        sensorEXP = '(</Sensor>)';
        if ~isempty(sensorStr)
            activeTF  = regexp(sensorStr, '(?<=Size=")\d(?=")', 'match');
            if  any(size(activeTF{1})==0) || char(activeTF{1})~='0'
                levl=1;
            end
            continue
        end
    else
% $$$         if regexp(tline,'</CalibrationCoefficients>')
% $$$             break
% $$$         end
        sensorEndStr = regexp(tline, sensorEXP);    
        if sensorEndStr
            levl=0;
            note=0;
        elseif levl==1
            theSensor  = regexp(tline, '(?<=<)[^\s]*', 'match');            
            theSensorID= regexp(tline, '(?<=SensorID=").*(?=")', 'match');
            sensorEXP  = ['(</',char(theSensor),'>)'];
            if  ismember(char(theSensor),fields(config)),
                disp('Duplicate sensor names, appending sensor ID to name')
                theSensor = {[theSensor{1},theSensorID{1}]};
            end
            configSensorSTR = ['config(1).',char(theSensor)];
            eval([configSensorSTR,'.ID=',char(theSensorID),';'])
            levl=2;
        elseif levl==2
            theToken = regexp(tline,'<(\w+)>.*</\1>','tokens');
            try theVar   = char(theToken{1});
                theValue = regexp(tline,['(?<=<',theVar,'>).*(?=</',theVar,'>)'],'match');
                if isempty(theValue), continue, end
                theValue = {['"',char(theValue{1}),'"']};
            catch theToken = regexp(tline,'<!--(.*)-->','tokens');
                try theValue = {['"',char(theToken{1}),'"']};
                theVar   = sprintf('note%d',note);
                note = note+1;
                catch theToken = regexp(tline,'(?<=Coefficients equation=")(\d)(?=")','match');
                    if coeff==0
                        configSensorSTRprefix = configSensorSTR;
                        configSensorSTR = [configSensorSTR, '.Coefficients',char(theToken)];
                        coeff=1;
                    elseif ~isempty(regexp(tline,'(</.*Coefficients>)'))
                        configSensorSTR = configSensorSTRprefix;
                        coeff=0;
                    end
                    continue
                end
            end
            eval([configSensorSTR,'.',theVar,'=',char(theValue{1}),';'])
        end
    end       
end

end