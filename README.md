# Hyperquicksort

## Example:

```
mpicc -O2 hyperquicksort.c -o hqs.o -lm
mpirun -np 4 ./hqs.o input_file res.csv -save
```
