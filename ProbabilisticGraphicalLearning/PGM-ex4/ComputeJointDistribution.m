function Joint = ComputeJointDistribution(F)

    numFactors = length(F);
    if numFactors<1
        Joint = struct('var',[],'card',[],'val',[]);
        return
    end
    
    Joint = F(1);
    for i=2:numFactors
        Joint = FactorProduct(Joint, F(i));
    end
    
end
