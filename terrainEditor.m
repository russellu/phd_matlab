clear all ; close all

cd C:\Users\butr2901\Pictures ; 
ground = imread('ground.png'); 

imagesc(ground); 
while true
   %[x,y] = ginput(1);  
   %disp(['x = ',num2str(x),' y=',num2str(y)]); 
   a = get(0, 'PointerLocation') ;
   %disp(a); 
    % check for keys
    k=get(gcf,'CurrentCharacter');
    if k~='@' % has it changed from the dummy character?
    set(gcf,'CurrentCharacter','@'); % reset the character
    % now process the key as required
    if k=='q', finish=true; end
    end
   pause(0.05); 
end