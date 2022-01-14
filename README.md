# Hyperquicksort

## Example:

To compile:
```
mpicc -O2 hyperquicksort.c -o hqs.o -lm
```
To execute:
```
mpirun -np 4 ./hqs.o input_file res.csv -save
```
