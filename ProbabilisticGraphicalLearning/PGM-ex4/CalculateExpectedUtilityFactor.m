% Copyright (C) Daphne Koller, Stanford University, 2012

function EUF = CalculateExpectedUtilityFactor( I )

  % Inputs: An influence diagram I with a single decision node and a single utility node.
  %         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
  %              the child variable = D.var(1)
  %         I.DecisionFactors = factor for the decision node.
  %         I.UtilityFactors = list of factors representing conditional utilities.
  % Return value: A factor over the scope of the decision rule D from I that
  % gives the conditional utility given each assignment for D.var
  %
  % Note - We assume I has a single decision node and utility node.
  EUF = [];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % YOUR CODE HERE...
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  R = I.RandomFactors;
  D = I.DecisionFactors;
  U = I.UtilityFactors;
  varD =[];
  for i=1:length(D)
    varD = [varD, D(i).var];
  end
  varD = unique(varD);

  %Joint = ComputeJointDistribution(F);
  %EUList = FactorProduct(Joint, U);
  %EU = sum(EUList.val); 
  %EUList =[];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %OptimizeWithJointUtility
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  JointUtilityFactor = ComputeSumDistribution(U);
  Joint = ComputeJointDistribution([R, JointUtilityFactor]);
  varJoint = Joint.var;
  [C, iA] = setdiff(varJoint, varD);
  varJoint = C;
  EUF = FactorMarginalization(Joint, varJoint);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %OptimizeLinearExpectations
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %EUList = struct([]);
  %for j = 1:length(U) 
  %  Joint = ComputeJointDistribution([R, U(j)]);
  %  varJoint = Joint.var;
  %  [C, iA] = setdiff(varJoint, varD);
  %  varJoint = C;
  %  EUList(j) = FactorMarginalization(Joint, varJoint);
  %end
  
  %EUF = ComputeSumDistribution(EUList);
  
end  
