cd C:\shared ;
vid = load('allframes.mat') ; 
vid = vid.allframes ; 

st = 764 ; en = size(vid,3)-2950 ; 
newvid = vid(60:626,31:374,st:en) ; 
mvid = uint8(mean(newvid,3)) ;
for i=1:size(newvid,3) ; newvid(:,:,i) = newvid(:,:,i)-mvid ; end
rotvid = zeros(size(newvid,2),size(newvid,1),size(newvid,3)) ; 
for i=1:size(newvid,3) ; rotvid(:,:,i) = imrotate(newvid(:,:,i),270) ; end
endwedge = 7095 ; 
