%% Calculating The Muo!
function MuoValue=CalculatingMuo(Datas,Rules,MFN,MFType,UpBnd,LowBnd)
switch MFType
    case 1              % For Triangular MF.
        Step = (UpBnd-LowBnd)/(MFN-1);
        Center = LowBnd+(Rules-1)*Step;
        MuoValue = trimf(Datas,[Center-Step,Center,Center+Step]);
        
    case 2              % For Trapezoid MF.
        Step = (UpBnd-LowBnd)/(MFN-1)/3;
        Center = LowBnd+(Rules-1)*3*Step;
        MuoValue = trapmf(Datas,[Center-2*Step,Center-Step,Center+Step,Center+2*Step]);        
        
    case 3              % For Gaussian MF.
        Step = (UpBnd-LowBnd)/(MFN-1)/2;
        Center = LowBnd+(Rules-1)*2*Step;
        MuoValue = gaussmf(Datas,[Step,Center]);
    otherwise
        disp(' Selected Membership Function is WRONG!')
        disp(' Please Start Again.')
end
end