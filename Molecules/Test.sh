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
cd ../C3H3NO
echo "############"
echo "10th test : not normal + Extreme Target"
echo "############"
../../source/DVCI C3H3NO_NotNormal.in
wait 
../../source/DVCI C3H3NO_5NU6.in 
