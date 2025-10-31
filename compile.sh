gcc main.c -o main -fopenmp -O3 \
    -Ibib -I./igraph/build/include -I./igraph/include \
    ./igraph/build/src/libigraph.a -llapack -lblas -lgfortran -lxml2 -lglpk -lgmp -lgsl -lgslcblas -lm -lstdc++

./main "./example/edges.txt" "./example/od.txt"
#./main "./fortaleza/edges_bairro/edges_fortaleza_0.txt" 50 0
#THREADS=8
#for i in $(seq 0 999); do
#    ./main "./fortaleza/edges_bairro/edges_fortaleza_0.txt" 50 $i &
#    while (( $(jobs -r | wc -l) >= THREADS )); do
#        sleep 1
#    done
#done