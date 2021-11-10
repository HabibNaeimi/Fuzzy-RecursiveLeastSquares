%% Creating the Rule Base.
function Rules=RuleBase(MFN,InputNum)

    PR = prod(MFN);
    Repeat = 1;
    Rules = zeros(PR,InputNum);
    for i=1:InputNum
        PR = PR/MFN(i);
        for l=1:Repeat
            
            for k=1:MFN(i)
                
                for j=1:PR
                    Rules(j+(k-1)*PR+(l-1)*PR*MFN(i),i) = k;
                end
                
            end
        end
        Repeat = Repeat*MFN(i);
    end
end