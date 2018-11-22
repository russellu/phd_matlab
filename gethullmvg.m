function hullmvg = gethullmvg(data,period)
halfperiod = round(period/2); 
wmvg1 = getwmvg(data,halfperiod)*2; 
wmvg2 = getwmvg(data,period); 
submvg = wmvg1-wmvg2; 
hullmvg = getwmvg(submvg,round(sqrt(period))); 
end

