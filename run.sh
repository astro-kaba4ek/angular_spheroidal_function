mkdir build 
cd build &&
cmake .. -DCMAKE_BUILD_TYPE=Release && 
make main && 
time ./main --input ../input.txt
