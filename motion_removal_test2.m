cd C:\shared\save_paper 

movie = load('movie_allmotions') ; movie = movie.allmotions ;

for i=1:8 ; subplot(2,4,i) ; imagesc(log(squeeze(movie(i,:,:)))) ; end ; 



for i=1:8

movie1 = (squeeze(movie(i,:,:))) ; binmovie = zeros(size(movie1)) ; 
for i=1:size(movie1,1)
    [z,mu,sigma] = zscore(movie1(i,:)) ; 
    binmovie(i,:) = movie1(i,:) > mu + sigma/2 ; 
end

figure,imagesc(binmovie) ; 

end