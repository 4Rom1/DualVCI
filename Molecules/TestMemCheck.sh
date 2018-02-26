#After running 
# $ ./TestMemCheck.sh > MemCheck.txt 2>&1
# The command 
# $ grep "ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)" MemCheck.txt | wc -l
# should return 11
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
cd ../C3H3NO
echo "############"
echo "9th Extreme Target"
echo "############"
valgrind --leak-check=full ../../source/DVCI C3H3NO/C3H3NO_5NU6.in 
