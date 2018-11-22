
function [OutputMat] = RusselFillOutside(InputMat)


row, column, z = size(InputMat);

OutputMat = zeros(size(InputMat));

for zIter = 1:z   
    for rowIter = 1:row
       
        rowInput = InputMat(rowIter, : , zIter);
        
        startId = find(diff(rowInput)==1);
        endId = find(diff(rowInput)==-1);
        
        newRow = zeros(size(rowInput));
        
        newRow(1:startId(1)) = 1;
        newRow(endId(end):end) = 1;
        
        OutputMat(rowIter, : , zIter) = newRow;
        
    end     
end