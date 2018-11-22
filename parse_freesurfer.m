function labs = parse_freesurfer 

subs = {'sub_alex','sub_charest','sub_esteban','sub_fabio','sub_gab','sub_gabriella','sub_genevieve','sub_gina','sub_jeremie','sub_julie','sub_katrine','sub_lisa'...
        ,'sub_marc','sub_marie','sub_mathieu','sub_maxime','sub_mingham','sub_patricia','sub_po','sub_russell','sub_sunachakan','sub_vincent'} ;
txtfiles = {
'lh_area.txt'
'rh_area.txt'
'lh_volume.txt'
'rh_volume.txt'
'lh_thickness.txt'
'rh_thickness.txt'
'lh_thicknessstd.txt'
'rh_thicknessstd.txt'
'lh_meancurv.txt'
'rh_meancurv.txt'
'lh_gauscurv.txt'
'rh_gauscurv.txt'
'lh_foldind.txt'
'rh_foldind.txt'
'lh_curvind.txt'
'rh_curvind.txt'
'aseg_volume.txt' 
'aseg_mean.txt'
'aseg_std.txt'} ;

for s=1:length(subs) ; 
    cd c:/shared/freesurfer ; sub = dir([subs{s},'-*']) ; 
    cd([sub.name,'/stats']) ; ls 
for txtfile=1:length(txtfiles) ;
    fid = fopen(txtfiles{txtfile}) ; 
    fline = 0 ; lcount = 1 ; 
    while fline ~= -1
       fline = fgets(fid) ;  
       disp(fline) ; 
       lines{lcount} = fline ; 
       l = lines{lcount} ; l(isspace(l)) = ' ' ; 
       splits{lcount} = strsplit(l,' ') ; 
       lsize = length(splits{1}) ; 
       if fline ~= -1
           if lcount==1 ; 
              labs{s,txtfile,lcount} = splits{lcount}(2:lsize-1) ; 
           elseif lcount==2
               labs{s,txtfile,lcount} = cellfun(@str2num,splits{lcount}(2:lsize-1)) ;
           end
       end
       lcount = lcount + 1 ; 
    end
    fclose(fid) ; 
end
end

end