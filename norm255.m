function normed = norm255(img)

normed = img - min(min(min(img))) ; 
normed = floor((normed./max(max(max(normed))))*255) ; 




end