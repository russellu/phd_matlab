function macd = getmacd(data, slowperiod, fastperiod)
% function macd = getmacd(data, slowperiod, fastperiod)

macd = getmvg(data,slowperiod) - getmvg(data,fastperiod);

end