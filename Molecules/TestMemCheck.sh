#After running 
# $ ./TestMemCheck.sh > MemCheck.txt 2>&1
# The command 
# $ grep "ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)" MemCheck.txt | wc -l
# should return 19
cd N2H2
echo "############"
echo "First test with wrong input file"
echo "############"
valgrind --leak-check=full ../../source/DVCI 
wait
valgrind --leak-check=full ../../source/DVCI Void.in
wait
valgrind --leak-check=full ../../source/DVCI N2H2_FF.in 
wait
echo "############"
echo "Second test with minimal input file"
echo "############"
valgrind --leak-check=full ../../source/DVCI N2H2_MIN.in 
wait
echo "############"
echo "Third test with different PESType"
echo "############"
valgrind --leak-check=full ../../source/DVCI N2H2_PESTYPE.in 
wait
echo "############"
echo "Fourth test with a target"
echo "############"
valgrind --leak-check=full ../../source/DVCI N2H2_Target.in 
wait
echo "############"
echo "Fifth test with KNREZ small"
echo "############"
valgrind --leak-check=full ../../source/DVCI N2H2_KNREZ_Small.in 
wait
echo "############"
echo "6th test with KNZREZ small"
echo "############"
valgrind --leak-check=full ../../source/DVCI N2H2_KNZREZ_Small.in 
wait
echo "############"
echo "7th test with KNNZ small"
echo "############"
valgrind --leak-check=full ../../source/DVCI N2H2_KNNZ_Small.in 
wait
echo "############"
echo "8th test with DoGraph=0"
echo "############"
valgrind --leak-check=full ../../source/DVCI N2H2_NoGraph.in 
wait
echo "############"
echo "9th test with ThrCoorSmall"
echo "############"
valgrind --leak-check=full ../../source/DVCI N2H2_ThrCoorSmall.in 
wait
echo "############"
echo "11th test with ThrCoorBig"
echo "############"
valgrind --leak-check=full ../../source/DVCI N2H2_ThrCoorBig.in 
wait
echo "############"
echo "12th test with target non reached"
echo "############"
valgrind --leak-check=full ../../source/DVCI N2H2TargetUnreached.in 
wait
echo "############"
echo "13th test with MaxEV big ThrCoor small"
echo "############"
valgrind --leak-check=full ../../source/DVCI N2H2_Target_MaxEV_Big.in 
wait
cd ../H2O
echo "############"
echo "14th test : Water"
echo "############"
valgrind --leak-check=full ../../source/DVCI H2O_Key.in
wait
echo "############"
echo "15th test : Water transitions"
echo "############"
#To check good behaviours, in Transitions.cc set IncK2 to 0, comment DoRot=0 after
#//DoRot must be set up to zero anyway, because Coriolis not part of operator
#Then the test will be done on the regular PES files. 
valgrind --leak-check=full ../../source/Transitions H2O_Key.in
wait #DVCI should be re-run before
valgrind --leak-check=full ../../source/DVCI H2O_AllIR.in
wait
valgrind --leak-check=full ../../source/Transitions H2O_AllIR.in
wait #DVCI should be re-run before
valgrind --leak-check=full ../../source/DVCI H2O_Target.in
wait
valgrind --leak-check=full ../../source/Transitions H2O_Target.in
