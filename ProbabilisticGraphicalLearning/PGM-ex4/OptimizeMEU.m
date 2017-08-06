% Copyright (C) Daphne Koller, Stanford University, 2012

function [MEU OptimalDecisionRule] = OptimizeMEU( I )

  % Inputs: An influence diagram I with a single decision node and a single utility node.
  %         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
  %              the child variable = D.var(1)
  %         I.DecisionFactors = factor for the decision node.
  %         I.UtilityFactors = list of factors representing conditional utilities.
  % Return value: the maximum expected utility of I and an optimal decision rule 
  % (represented again as a factor) that yields that expected utility.
  
  % We assume I has a single decision node.
  % You may assume that there is a unique optimal decision.
  %D = I.DecisionFactors(1);
  D = I.DecisionFactors;
  U = I.UtilityFactors;
  R = I.RandomFactors;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % YOUR CODE HERE...
  % 
  % Some other information that might be useful for some implementations
  % (note that there are multiple ways to implement this):
  % 1.  It is probably easiest to think of two cases - D has parents and D 
  %     has no parents.
  % 2.  You may find the Matlab/Octave function setdiff useful.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  EUF = CalculateExpectedUtilityFactor(I);
  DecisionRules = ComputeJointDistribution(D);
  numD = length(DecisionRules.var);
  varList = zeros(1,numD);
  cardList = zeros(1,numD);
  for i=1:numD
    varList(i) = DecisionRules.var(i);
    cardList(i) =DecisionRules.card(i);
  end
  [varD, map]=unique(varList);
  cardD = cardList(map); 
  
  varChD = zeros(1,length(D));
  varPaD = [];

  for i=1:length(D)
    varChD(i) = D(i).var(1);
  end

  varChD = unique(varChD);
  varPaD = unique(setdiff(varD, varChD));
  
  %mapCh = arrayfun(@(x) find(varD==x),varChD);
  %mapPa = arrayfun(@(x) find(varD==x),varPaD);
  [dummy, mapCh] = ismember(varChD,varD);
  [dummy, mapPa] = ismember(varPaD,varD);
  cardChD = cardD(mapCh);
  cardPaD = cardD(mapPa);
  
  %[MEU,idx] = max(EUF.val);
  %assign = IndexToAssignment(idx,cardD);
  %assign(mapCh)
  OptimalDecisionRule = DecisionRules;
  OptimalDecisionRule.val = zeros(1,length(DecisionRules.val));
  assign = zeros(1,numD);
  utility = zeros(prod(cardChD),prod(cardPaD));
  %maxUtility = zeros(prod(cardPaD));
  %maxJ = zeros(prod(cardPaD));
  MEU = 0;
  for i=1:prod(cardPaD)
    assignPa = IndexToAssignment(i,cardPaD);
    assign(mapPa) = assignPa;
    for j=1:prod(cardChD)
      assignCh = IndexToAssignment(j,cardChD);
      assign(mapCh) = assignCh;
      utility(j,i) = EUF.val(AssignmentToIndex(assign,cardD));
    end
    [maxUtility,maxJ] = max(utility(:,i));
    MEU = MEU + maxUtility;
    maxCh = IndexToAssignment(maxJ,cardChD);
    assign(mapCh) = maxCh;
    idx = AssignmentToIndex(assign,cardD);
    OptimalDecisionRule.val(idx) = 1;
  end
  
      
      
      
end
