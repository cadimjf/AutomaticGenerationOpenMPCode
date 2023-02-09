g++ -I/usr/local/include/ -I/usr/local/include/cvode -I/usr/local/include/sundials/ -I/usr/local/include/nvector -L /usr/local/lib -o main main.cpp Stopwatch.cpp -lsundials_cvode -lsundials_nvecserial -lm -w -O3 -fopenmp

# echo "----------------"
./main 1
echo "-----------------"
# ./main 2
# echo "-----------------"
# ./main 3
# echo "-----------------"
# ./main 4
# echo "-----------------"
# ./main 5
# echo "-----------------"
