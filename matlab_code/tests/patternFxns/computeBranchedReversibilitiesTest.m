classdef computeBranchedReversibilitiesTest < matlab.unittest.TestCase
    
    methods (Test)
        function testComputeBranchedReversibilities1(testCase)
            
            revMatrix = [0,1,0,1,1,0,1,0,1;0,1,0,1,1,1,0,1,0;1,0,1,0,1,0,1,0,1;1,0,1,0,1,1,0,1,0];
            lbRev = [0.00400964663161663;0.151021586158862;0.0894718942001685;...
                     0.0189704801904070;0.0865906933913647;0.0753094930336852;...
                     0.133242559631399;0.184404535137972;0.256979111624525];
                 
            branchedRev = computeBranchedReversibilities(revMatrix,lbRev);
            
            trueRes = [0.00400964663161663;0.151021586158862;0.165982419717652;...
                       0.0189704801904070;0.0865906933913647;0.0753094930336852;...
                       0.133242559631399;0.668107747225681;0.610174680627967];
           
            testCase.verifyLessThanOrEqual(abs(trueRes - branchedRev), abs(trueRes * 10^-14));
            
        end
     
    end
end

