function [inspiredV,expiredV] = BreathIntegrator(inlocs,exlocs,flowScaledFiltered)

%get the integral of the first few cycles of data
%if integral > .10
%shift the baseline slightly positive
%else if integral < 0.1
%shift the baseline slightly negative
%else done.

%use locs to get the index of each peak
%from locs decrement the index until mag is < 0 record index - lower,
%increment the index until mag < 0 record index - higher
%integrate betwen lower and upper index
 inspiredV = zeros(length(inlocs),1);
 for i = 1:length(inlocs)
  lower = inlocs(i);
  while (flowScaledFiltered(lower) > 0)
   lower = lower - 1;
  end
  lower = lower + 1;
  upper = inlocs(i);
  while (flowScaledFiltered(upper) > 0)
   upper = upper + 1;
  end
  upper = upper - 1;
  inspiredV(i) = trapz(0.01,flowScaledFiltered(lower:upper));
 end
 inspiredVT = transpose(inspiredV);

 expiredV = zeros(length(exlocs),1);
 for i = 1:length(exlocs)
  lower = exlocs(i);
  while (flowScaledFiltered(lower) < 0)
   lower = lower - 1;
  end
  lower = lower + 1;
  upper = exlocs(i);
  while (flowScaledFiltered(upper) < 0)
   upper = upper + 1;
  end
  upper = upper - 1;
  expiredV(i) = trapz(0.01,flowScaledFiltered(lower:upper));
 end 
 expiredVT = transpose(expiredV);
 
 flowScaledFilteredmean = mean(flowScaledFiltered(inlocs(1):exlocs(end)));
 inexVdiff = mean(inspiredV) + mean(expiredV);
 flowScaledFiltered = flowScaledFiltered - flowScaledFilteredmean; 
end