
To obtain the V_ks matrix
in our case, V_k is a column in V_ks such that its columns are indexed by strike_id

>>generate_hyperparameters(100,(44:1:56)',0.1,0.01,0,0.25,0,2^6)
>>load("hyperparamters")
>>V_ks = IV_test_v1(S_0,strikes,T,r,q,c1,c2,c4,numGrid,1,isHes)