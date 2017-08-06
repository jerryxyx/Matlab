function Sum = ComputeSumDistribution(F)

    numFactors = length(F);
    if numFactors<1
        Sum = struct('var',[],'card',[],'val',[]);
        return
    end
    
    Sum = F(1);
    for i=2:numFactors
        Sum = FactorAdd(Sum, F(i));
    end
    
end
