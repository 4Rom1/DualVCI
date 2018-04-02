cd N2H2
echo "############"
echo "First test with wrong input files"
echo "############"
../../source/DVCI 
wait
../../source/DVCI Void.in
wait
../../source/DVCI N2H2_FF.in 
wait
../../source/DVCI N2H2_NotNormal.in 
wait
echo "############"
echo "Second test with minimal input file"
echo "############"
../../source/DVCI N2H2_MIN.in 
wait
echo "############"
echo "Third test with different PESType"
echo "############"
../../source/DVCI N2H2_PESTYPE.in 
wait
echo "############"
echo "Fourth test with a target"
echo "############"
../../source/DVCI N2H2_Target.in 
wait
echo "############"
echo "Fifth test with KNREZ small"
echo "############"
../../source/DVCI N2H2_KNREZ_Small.in 
wait
echo "############"
echo "6th test with KNZREZ small"
echo "############"
../../source/DVCI N2H2_KNZREZ_Small.in 
wait
echo "############"
echo "7th test with KNNZ small"
echo "############"
../../source/DVCI N2H2_KNNZ_Small.in 
wait
echo "############"
echo "8th test with DoGraph=0"
echo "############"
../../source/DVCI N2H2_NoGraph.in 
wait
echo "############"
echo "9th test with DeltaNEV"
echo "############"
../../source/DVCI N2H2_DeltaNEV.in 
wait
echo "############"
echo "10th test with ThrCoorSmall"
echo "############"
../../source/DVCI N2H2_ThrCoorSmall.in 
wait
echo "############"
echo "11th test with ThrCoorBig"
echo "############"
../../source/DVCI N2H2_ThrCoorBig.in 
wait
echo "############"
echo "12th test with target non reached"
echo "############"
../../source/DVCI N2H2TargetUnreached.in 
wait
echo "############"
echo "13th test with MaxEV big ThrCoor small"
echo "############"
../../source/DVCI N2H2_Target_MaxEV_Big.in 
wait
echo "############"
echo "14th test with Kappa Small"
echo "############"
../../source/DVCI N2H2_KappaSmall.in
wait
echo "############"
echo "15th test with MaxEV=0"
echo "############"
../../source/DVCI N2H2_MAXEV0.in
wait
echo "############"
echo "16th test with MaxFreq=3000"
echo "############"
../../source/DVCI N2H2_3000.in
wait
cd ../C3H3NO
echo "############"
echo "17th test : not normal + Extreme Target"
echo "############"
../../source/DVCI C3H3NO_NotNormal.in
wait 
../../source/DVCI C3H3NO_5NU6.in 
wait
cd ../H2O
echo "############"
echo "18th test : Water"
echo "############"
#Fail on Ground state.
../../source/DVCI H2O_FailG.in
#Simple fundamentals calculation.
../../source/DVCI H2O_Key.in
wait
echo "############"
echo "19th test : Water transitions"
echo "############"
#To check good behaviours, in Transitions.cc set IncK2 to 0, comment DoRot=0 after
#//DoRot must be set up to zero anyway, because Coriolis not part of operator
#Then the test will be done on the regular PES files. 
../../source/Transitions H2O_Key.in
wait #DVCI should be re-run before
../../source/DVCI H2O_AllIR.in
wait
../../source/Transitions H2O_AllIR.in
wait #DVCI should be re-run before
../../source/DVCI H2O_Target.in
wait
../../source/Transitions H2O_Target.in
wait
../../source/DVCI H2O_Fail1.in
wait
../../source/Transitions H2O_Fail1.in
echo "############"
echo "20th test : Water FinalVCI"
echo "############"
../../source/DVCI H2O_Fund.in
wait
../../source/FinalVCI H2O_Targ.in
wait #DVCI should be re-run before
../../source/DVCI H2O_AllIR2.in
wait
../../source/FinalVCI H2O_AllIR2.in
wait #DVCI should be re-run before
../../source/DVCI H2O_Target2.in
wait
../../source/FinalVCI H2O_Target2.in
wait
../../source/DVCI H2O_Fail2.in
wait #Fail because not exit success for DVCI
../../source/FinalVCI H2O_Fail2.in
wait #No Rotational terms
../../source/DVCI H2O_AllIRNoRot.in
wait
../../source/FinalVCI H2O_AllIRNoRot.in
