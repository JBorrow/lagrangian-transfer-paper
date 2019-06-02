for radius in 1.0 1.00001 1.1 1.2 1.3 1.4 1.5;
do
    cd $radius

    mkdir plots
    cd plots

    cp ../../plots/* .
    sbatch  --dependency=afterok:941994 submit.slurm

    cd ../..
done

