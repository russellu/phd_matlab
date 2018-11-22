clear all ;close all ; 
subs = {'alex','charest','esteban','fabio','gab','gabriella','genevieve','gina','guillaume','jeremie','julie','katrine','lisa','marc',...
    'marie','mathieu','maxime','mingham','patricia','po','russell','sunachakan','tah','vincent'} ; 
stimnames = {'unperturbed','5% contrast','33% contrast','plaid','10% random','60% random'};

for sb=1:length(subs); disp(sb) ; 
    cd(['e:/allfmris/sub_',subs{sb}]) ;  
    clear allstims; 
    commons = dir('common*gz');
    corrs = load_untouch_nii('meancorrs.nii.gz') ; 
    
    corrs.img(:,1:40,:) = 0 ;

    [sv,si] = sort(corrs.img(:),'ascend') ;
    
    stimcounts = zeros(1,6);
    voxels = si(1:100) ;
    [vx,vy,vz] = ind2sub(size(corrs.img),voxels) ; 
    for cm=1:length(commons)
        nii = load_untouch_nii(commons(cm).name);
        
        for i=1:length(vx)
           ts(i,:) = squeeze(nii.img(vx(i),vy(i),vz(i),:));   
        end
        %stdts = std(zscore(ts,[],2),0,2); goods = find(zscore(stdts)<1) ; 
        meants = mean(ts(:,:)) ;
        cmt = load(['trigs/stimTimes',num2str(cm)]) ; cmt = cmt.stimTimes; cmt = cell2mat(cmt); 
        stypes = cmt(2:2:end) ; stimes = cmt(1:2:end) ; 
        strs = round(stimes/2) ; 
        for i=1:length(stypes)
            allstims(stypes(i),stimcounts(stypes(i))+1,:) = meants(strs(i)-2:strs(i)+8) ;
            stimcounts(stypes(i)) = stimcounts(stypes(i))+1;
        end
    end
    percstims(sb,:,:,:) = (allstims - repmat(mean(allstims(:,:,3),3),[1,1,11])) ./ repmat(mean(allstims(:,:,3),3),[1,1,11]) ;
    neg_mstims = squeeze(mean(percstims(sb,:,:,:),3));
    save('neg_mstims','neg_mstims');         
end

m_neg = squeeze(mean(percstims,3)); 
save('m_neg','m_neg'); 
