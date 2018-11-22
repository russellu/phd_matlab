function pstr = format_p(pval)
%function pstr = format_p(pval)

if pval < 0.0001
    pstr = 'p<0.0001';  
elseif pval < 0.001
    pstr = 'p<0.001';     
elseif pval < 0.01
    pstr = ['p=',num2str(round(pval*1000)/1000)]; 
elseif pval < 0.05
    pstr = ['p=',num2str(round(pval*1000)/1000)];   
else
    pstr = ['p=',num2str(round(pval*1000)/1000)];  
end



end