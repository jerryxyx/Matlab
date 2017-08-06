% Copyright (C) Daphne Koller, Stanford University, 2012

function [MEU OptimalDecisionRule] = OptimizeLinearExpectations( I )
  % Inputs: An influence diagram I with a single decision node and one or more utility nodes.
  %         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
  %              the child variable = D.var(1)
  %         I.DecisionFactors = factor for the decision node.
  %         I.UtilityFactors = list of factors representing conditional utilities.
  % Return value: the maximum expected utility of I and an optimal decision rule 
  % (represented again as a factor) that yields that expected utility.
  % You may assume that there is a unique optimal decision.
  %
  % This is similar to OptimizeMEU except that we will have to account for
  % multiple utility factors.  We will do this by calculating the expected
  % utility factors and combining them, then optimizing with respect to that
  % combined expected utility factor.  
  MEU = [];
  OptimalDecisionRule = [];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % YOUR CODE HERE
  %
  % A decision rule for D assigns, for each joint assignment to D's parents, 
  % probability 1 to the best option from the EUF for that joint assignment 
  % to D's parents, and 0 otherwise.  Note that when D has no parents, it is
  % a degenerate case we can handle separately for convenience.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  D = I.DecisionFactors;
  U = I.UtilityFactors;
  R = I.RandomFactors;
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
