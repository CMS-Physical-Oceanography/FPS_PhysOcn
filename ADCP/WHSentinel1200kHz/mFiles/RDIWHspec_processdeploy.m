 function [Data] = RDIWHspec_processdeploy(path)

   %======================================================================
   % This function inputs the path of a deployment folder from Wavemon
   % processing. The spectrum files stored in the folder are read into the
   % file. It outputs a 64xNbursts array of forier coeffs
   % from each burst in the deployment, and the timestamp of each burst. 
   %======================================================================
     % read contents of path
       source={fullfile(path,'PSpec*'),fullfile(path,'SSpec*'), fullfile(path,'VSpec*')};
       Data=struct();
       bur_names = ['P_burst';'S_burst';'V_burst'];
       time = ['P_time'; 'S_time';'V_time'];
       f = ['P_freq'; 'S_freq';'V_freq'];
       for i=1:length(source)
           files=dir(source{i});
           files={files.name};
           % initialize output arrays
           dim = size(files); 
           % initial timestamp array
           times = strings(1,dim(2));
           % initial burst array 
           bursts = zeros(64,dim(2));
                               
        for j=1:dim(2)
           file = files{1,j};
           fil=append(path,'\', file);
           strdata = fileread(fil);
           % slice data from header for Vspec, Pspec, Spec files 
           data0idx = strfind(strdata,'830078)'); 
           data = str2num(strdata(data0idx+14:end));
           burst = data;
           bursts(:,j)=burst;
           num_expression = '(\d+)';
           match(j) = regexp(file, num_expression, 'match');
           times(j)=datetime(match{j},'InputFormat','yyMMddHHmm');
         
        end
          %Reading files for frequency bands and units
           fid=fopen(fil); 
           dum=fgetl(fid);
           tline=fgetl(fid);
           units = regexp(tline,'(?<=units of\s).*','match');
          %Andvance to the third line of the .txt files to  read freq bands
           tline=fgetl(fid);
           f_incr = str2double(regexp(tline,'(?<=are\s).*(?=\sHz)','match'));
           f_unit=regexp(tline,'(?<=0.01562500\s).*(?=\swide)','match');
           f0= str2double(regexp(tline,'(?<=at\s).*(?=\))','match'));
           num_freq= str2double(regexp(tline,'(?<=%\s).*(?=\sFrequency)','match'));
           freq=f0:f_incr:(num_freq)*f_incr;
       %Savings burst and time to a variable before being overwritten by
       %the next loop.
        eval(sprintf('%s = bursts;',bur_names(i,:)))
        eval(sprintf('%s = times;',time(i,:)))
        eval(sprintf('%s = freq;',f(i,:)))
        clear bursts; clear times; clear freq;
        disp('done')
       end
       
       
       %DSpec File
       source=fullfile(path, 'DSpec*');
        files=dir(source);
           files={files.name};
           fil=append(path,'\', files);
           % initialize output arrays
           dim = size(files); 
           % initial timestamp array
           times = strings(1,dim(2));
           % initial burst array 
           bursts = zeros(64,90,dim(2));
          %Reading files for spectrum units
           fid=fopen(fil{1}); 
           dum=fgetl(fid); dum=fgetl(fid); dum=fgetl(fid);
           tline=fgetl(fid);
           dir_units = regexp(tline,'(?<=are\s).*(?=\per cycle)','match');
                  
        for k=1:dim(2)
           file = files{1,k};
           fil=append(path,'\', file);
           strdata = fileread(fil);
           dir0idx = strfind(strdata,'begins at'); % find index 
           dir0 = strdata(dir0idx+10:dir0idx+12); % slice first bin edge
           dbinedges= str2num(dir0)+[0:89]*4;
           dbinedges(dbinedges>360)= dbinedges(dbinedges>360)-360;
           
           % slice data from header 
           data0idx = strfind(strdata,'degrees'); 
           data = str2num(strdata(data0idx+10:end));
           burst = data;
           bursts(:,:,i) = burst;
           %converting file name to datetime expression
           num_expression = '(\d+)';
           match(j) = regexp(file, num_expression, 'match');
           times(j)=datetime(match{j},'InputFormat','yyMMddHHmm');
          %Savings burst and time to a variable.
           eval(sprintf('%s = bursts;', 'D_burst'));
           eval(sprintf('%s = times;','D_time'));
        end
        
     %Data Structure for the four different spectrums 
     Data.PSpec.burst=P_burst;
     Data.PSpec.units=units;
     Data.SSpec.burst=S_burst;
     Data.SSpec.units=units;
     Data.VSpec.burst=V_burst;
     Data.VSpec.units=units;
     Data.DSpec.burst=D_burst;
     Data.DSpec.units=dir_units;
     Data.SSpec.dir_bins=dbinedges;
     Data.time=times;
           
 end
            
        
        
        
               
               