
%create a list of all file names...
filelist= dir ('C:\Users\wanho\Documents\Honours_Project\Duo_Daq_Data')
filelist(1).name

%have ti loop through this one
cd 'C:\Users\wanho\Documents\Honours_Project\R_subserialdata'

A= readtable('temp_humidity_data_R_1.csv');
Tmeans=[]
RHmeans=[]
Tmeans= A(:,2);
RHmeans= A(:,4);
Tmeans = table2array(Tmeans)
RHmeans = table2array(RHmeans)


RHmeans

cd 'C:\Users\wanho\Documents\Honours_Project\Duo_Daq_Data'

for i = 1:47
    
    cd 'C:\Users\wanho\Documents\Honours_Project\Duo_Daq_Data'

    if (i==1) || (i==2)
        continue

    end

    filename=filelist(i).name


    G = readtable(filename);  

    timestamps = G(:,1)
    pressure = G(:,2)
    Flow=G(:,3)
    
    rawpressurearray = table2array(pressure)
    rawflowarray = table2array(Flow)
    %define butterworth and apply
    [b,a] = butter(10,5/50,'low');
    filteredpressurearray = filter(b,a,rawpressurearray);


    %on the first go, the data is nowhere near the set ipap and epap
    %values, then i tried to see if the offset is connsistent and it turns
    %out when I divide it by 1.6 it resembles the set ipap and epap values.
    filteredpressurearray= filteredpressurearray/1.6
    
    [inpkspressure,inlocs] = findpeaks(filteredpressurearray,'MinPeakDistance',5,'MinPeakProminence', 0.2*max(filteredpressurearray(100:end)));
    
    invertedPressure = -filteredpressurearray
    [expks,exlocs] = findpeaks(invertedPressure,'MinPeakDistance',5,'MinPeakProminence', 0.2*max(filteredpressurearray(100:end)));
    
    %get the  last 20 IPAP
    trimmedIPAP = inpkspressure(end-20:end-1);
    inlocs= inlocs(end-20:end-1)
    
    %desciptive stats 
    IPAPmean = mean(trimmedIPAP);
    IPAPsd = std(trimmedIPAP);
    IPAPn = length(trimmedIPAP);
    
    %get the last 20 EPAP
    trimmedEPAP = -expks(end-20:end-1);
    exlocs= exlocs(end-20:end-1)
    
    %desciptive stats 
    EPAPmean = mean(trimmedEPAP);
    EPAPsd = std(trimmedEPAP);
    EPAPn = length(trimmedEPAP);
    
    
    
    
     
    
    masternpoint=[]
    peakpoints= inlocs
    
    
     %finding pressure/time product
     %I locate the index where it is 10 percent of (IPAPmean -EPAP mean),
     %
     % trough index .........
     %then this index plus 30 will be the index for 300ms and +50 will be index of
     %500ms.
     threehundredpressuretime=[]
     fivehundredpressuretime=[]
     for i = 1:length(peakpoints)
         
        
    
         indexpeak= peakpoints(i)
         requiredninetyn=0
         requiredtenn=0
         
         while indexpeak>0
             
             display (indexpeak)
             display(filteredpressurearray(indexpeak))

             
             % this is irrelavant 
             if (filteredpressurearray(indexpeak)- (0.9 * (IPAPmean -EPAPmean) + EPAPmean)  <0) && (requiredninetyn==0)
                 requiredninetyn= indexpeak
            
             end


             if (filteredpressurearray(indexpeak)- (0.1 * (IPAPmean - EPAPmean) + EPAPmean) <0)
                 requiredtenn= indexpeak+1
                 
                 %irrelavant
                 diffnpoints= requiredninetyn-requiredtenn
                 %irrelevant
                 masternpoint=[masternpoint,diffnpoints];
                 %irrelevant
                 risetimeofeachbreath= masternpoint*0.01
                 
                 %300ms = 30 interval 
                 %500ms = 50 interval 
                 threehundredpressuretime(end+1)=trapz(0.01,filteredpressurearray(requiredtenn:requiredtenn+30))
                 fivehundredpressuretime(end+1)=trapz(0.01,filteredpressurearray(requiredtenn:requiredtenn+50))
    
                 
    
    
    
                 break
       
    
             end
    
             
             
                
    
             
    
             indexpeak= indexpeak-1
    
    
         end
         i=i+1
          
    
     end 
    
     combinedpressuretime= [threehundredpressuretime;fivehundredpressuretime]
     combinedpressuretime= transpose(combinedpressuretime)
    
     %flow calculations:
     
     fflow = filter(b,a,rawflowarray);
     [pks,locs] = findpeaks(fflow,'MinPeakDistance',5,'MinPeakProminence', 0.9*max(fflow(100:end)));
     invertedFFlow = -fflow;
     [expks,exlocs] = findpeaks(invertedFFlow,'MinPeakDistance',5,'MinPeakProminence', 0.9*max(invertedFFlow(100:end)));
    

    
    


    
    %convert flow to voltage:
    %positive
    %Slope:156.00000
    %Intercept:-0.32800
    
    %negative 
    %Slope:154.00000
    %Intercept:0.43200
    offset1=-0.32800
    gain1=156.00000
    
    offset2=0.43200
    gain2=154.00000
    flowVolt = (rawflowarray - offset1) / gain1
    fplus = flowVolt > 0;%find index of all values greater than zero in flowVolt array
    flowScaled = zeros(length(rawflowarray),1);
    flowScaled(fplus) = gain1 * flowVolt(fplus) + offset1;%scale positive inspiratory values 
    flowScaled(~fplus) = gain2 * flowVolt(~fplus) + offset2;%scale negative expiratory values
    
    %filter flow data. InFlowd IPAP max and EPAP mins. 
    fflowVolt = filter(b,a,flowVolt);
    fflowVoltMean = mean(fflowVolt);
    flowVoltMean = mean(flowVolt);
    
    
    %filter scaled flow and find max inspiratory flow mean,sd and max expiratory flow mean,sd  
    flowScaledFiltered = filter(b,a,flowScaled);
    
    flowScaledFilteredMean = mean(flowScaledFiltered);
    flowScaledFiltered = flowScaledFiltered - flowScaledFilteredMean;
    
    [inpks,inlocsflow] = findpeaks(flowScaledFiltered,'MinPeakDistance',5,'MinPeakProminence', 0.6*max(flowScaledFiltered(100:end)));
    invertedflowScaledFiltered = -flowScaledFiltered;
    [expksflow,exlocsfinal] = findpeaks(invertedflowScaledFiltered,'MinPeakDistance',5,'MinPeakProminence', 0.6*max(invertedflowScaledFiltered(100:end)));

    

    inpks = inpks(end-20:end-1);
    inlocsflow = inlocsflow(end-20:end-1);

    expksflow = expksflow(end-20:end-1);
    exlocsfinal = exlocsfinal(end-20:end-1);
    exlocsfinal=transpose(exlocsfinal)
    
    
    %calculate volume by identifying peak and then go left and right until 0
    %to get upper and lower bound then integrate, i just copied and paste 

   
    
     
    
     %Inspiratory Calibration Factor = 0.0189
    %Expiratory Calibration Factor = -0.0186
    
    
    
    
    % for this part, I keep running into errors where you calculate the
    % difference between above zero and below zero and shift the baseline
    % 
    

     diffVOK = 0;
     diffVOKIt = 0;
     baselineAdj = 0;
     while (diffVOK ~= 1)
      cd 'C:\Users\wanho\Documents\Honours_Project'
      [volIn,volEx] = BreathIntegrator(inlocsflow,exlocsfinal,flowScaledFiltered);
      volIn = volIn * 0.0189;
      volEx = volEx * -0.0186;
      volInMean = mean(volIn);
      volExMean = mean(volEx);
      volMeanDiff = volInMean - volExMean;
      if (abs(volMeanDiff) > 0.003)
       baselineAdj = volMeanDiff * 0.05;
      %elseif (volMeanDiff < -0.003)
      % baselineAdj = volMeanDiff * -0.05;
      else
       diffVOK = 1;
      end
      flowScaledFiltered = flowScaledFiltered - baselineAdj;
      diffVOKIt = diffVOKIt + 1;
     end
     %&&&&&&&&&&&&& INSERT VOL IN CODE HERE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     
     %&&&&&&&&&&&&& INSERT VOL IN CODE HERE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     volInMean = mean(volIn * 1000);
     volInN = length(volIn);
     volExMean = mean(volEx * 1000);
     volExN = length(volEx);
     volInSD = std(volIn);
     volExSD = std(volEx);


   
    flowinmax= inpks
    flowinmaxsd= std(flowinmax)
    flowinmaxn= length(flowinmax)
    flowoutmax= expksflow
    flowoutmaxsd= std(flowoutmax)
    flowoutmaxn= length(flowoutmax)


    %%humidity stuff
    atmP = 766.0;%based on BOM MSLP for data for 24-11-2020 mean of 09:00 and 15:00 readings
    atmT = 25;
    atmH = 48;
    Pavs = 0.61121 * exp((18.678 - atmT/234.5) * (atmT / (257.14 + atmT))) * 7.50062;%Pav = Vapour pressure of H2O in ambient air
     %STPD = VATPS x (273 / (273 + T°C) x (Patm – Pdh2o)/760);Patm = 760mmHg
     Pah2os = atmH / 100 * Pavs;
     Vvstpd = volIn * (273 / (273 + atmT) * (atmP - Pah2os) / 760);
     
     Pdvs = 0.61121 * exp((18.678 - Tmeans(i-2) / 234.5) * (Tmeans(i-2) / (257.14 + Tmeans(i-2)))) * 7.50062;
     Pdh2os = RHmeans(i-2) / 100 * Pdvs;
     Vdstpd = volInMean * (273 / (273 + Tmeans(i-2)) * (760 - Pdh2os) / 760);
    
    

    M = zeros(24,length(trimmedIPAP));
    
    M(1,1:length(trimmedIPAP)) = trimmedIPAP;
    M(2,1:length(IPAPmean)) = IPAPmean;
    M(3,1:length(IPAPsd)) = IPAPsd;
    M(4,1:length(IPAPn)) = IPAPn;
    
    
    M(5,1:length(trimmedEPAP)) = trimmedEPAP;
    M(6,1:length(EPAPmean)) = EPAPmean;
    M(7,1:length(EPAPsd)) = EPAPsd;
    M(8,1:length(EPAPn)) = EPAPn;
    
    M(9,1:length(threehundredpressuretime)) = threehundredpressuretime;
    M(10,1:length(fivehundredpressuretime)) = fivehundredpressuretime;
    
    M(11,1:length(flowinmax)) = flowinmax
    M(12,1:length(flowinmaxsd)) = flowinmaxsd;
    M(13,1:length(flowinmaxn)) = flowinmaxn;
    
    M(14,1:length(flowoutmax)) = flowoutmax
    M(15,1:length(flowoutmaxsd)) = flowoutmaxsd;
    M(16,1:length(flowoutmaxn)) = flowoutmaxn;
    
    
    M(17,1:length(volIn)) = volIn;
    M(18,1:length(volInMean)) = volInMean;
    M(19,1:length(volInsd)) = volInsd;
    
    M(20,1:length(volEx)) = volEx;
    M(21,1:length(volExMean)) = volExMean;
    M(22,1:length(volExsd)) = volExsd;
    M(23,1:length(Vvstpd)) = Vvstpd;
    M(24,1:length(Vdstpd)) = Vdstpd;
    
    M = transpose(M)
    
    T = array2table( M ) 
    
    T.Properties.VariableNames = {'Trimmed IPAP','IPAP mean','IPAP SD ','IPAPN','Trimmed EPAP','EPAP mean','EPAP SD','EPAP n ','300ms risetime','500ms risetime','flowinmax','flow in maxsd','flow in max n ','flowout max','flow out max sd ','flow out max n','inspired vt','inspired vt sd','invol mean','expiredvt','exvol mean','ex vol sd','Vvstpd','Vdstpd'};
    
    cd 'C:\Users\wanho\Documents\Honours_Project\matlab_duo_daq_data'
    writetable(T,filename,'Delimiter',',')  

       
end
