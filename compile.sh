gcc main.c -o main -fopenmp -O3 \
    -Ibib -I./igraph/build/include -I./igraph/include \
    ./igraph/build/src/libigraph.a -llapack -lblas -lgfortran -lxml2 -lglpk -lgmp -lgsl -lgslcblas -lm -lstdc++

./main "./example/edges.txt" "./example/od.txt" "Dial"