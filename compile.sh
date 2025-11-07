gcc main.c -o main -fopenmp -O3 \
    -Ibib -I./igraph/build/include -I./igraph/include \
    ./igraph/build/src/libigraph.a \
    -llapack -lblas -larpack -lgfortran -lxml2 -lglpk -lgmp -lgsl -lgslcblas -lm -lstdc++


./main "./fortaleza/edges.txt" "iTAPAS" 10 300