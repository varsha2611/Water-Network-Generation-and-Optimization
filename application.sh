#PBS -N gnu-parallel-example
#PBS -qsub -I -l  select=1:ncpus=20:mem=120gb,walltime=72:00:00

module add gnu-parallel

cd $PBS_O_WORKDIR

Input_file = $1
Params = $2
python generate_inputs.py -i $1 -p $2

module add gnu-parallel
parallel -j100 <inputs.txt
