
rm -f strain_stress.txt
rm -f test_concrete.out

make 

./test_concrete.out

python plot_sol.py

