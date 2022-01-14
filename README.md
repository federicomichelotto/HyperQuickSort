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

Where:
- *input_file*: input file path.
- 'res.csv': csv file path in which to store the sorting time.
- '-save': flag to save the sorted list as a file 'input_file_sorted'.

