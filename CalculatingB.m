%% Calculating The B!
function FinalB=CalculatingB(Data,Rules,MFNum,MFType,UpBnd,LowBnd)

MuoValue = zeros(size(Rules));

    for i=1:size(Rules,1)
        for j=1:size(Rules,2)
            MuoValue(i,j) = CalculatingMuo(Data(j),Rules(i,j),MFNum(j),MFType(j),UpBnd(j),LowBnd(j));
        end
    end
    
    a = prod(MuoValue,2);
    b = sum(a);
    FinalB = a/b;
end