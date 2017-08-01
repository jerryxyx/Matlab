function testInference

  load PA3Models.mat;
  load PA3Data.mat;
  %printf("torturing");
  
  imageModel.ignoreSimilarity=true;
  wordPrediction = RunInference(BuildOCRNetwork(allWords{1}, imageModel, [], []))%rorruring  
  wordPrediction = RunInference(BuildOCRNetwork(allWords{1}, imageModel, pairwiseModel, tripletList))%rorturing
  ScoreModel({allWords{11:30}},imageModel,[],[])


  imageModel.ignoreSimilarity=false;
  wordPrediction = RunInference(BuildOCRNetwork(allWords{1}, imageModel, [], []))%rorruring exactly the same as the model ignoring similarity
  wordPrediction = RunInference(BuildOCRNetwork(allWords{1}, imageModel, pairwiseModel, tripletList))%torturing
  ScoreModel({allWords{11:30}},imageModel,[],[])%exactly the same as the model ignoring similarity
end