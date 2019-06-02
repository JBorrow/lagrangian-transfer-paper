for radius in 1.0 1.00001 1.1 1.2 1.3 1.4 1.5;
do
    mkdir $radius
    cd $radius

    cp ../lt/submit.slurm .
    cp ../lt/parameters.yml .

    cp ../lt/*.py .

    sed -i "s/RADIUS_PLACEHOLDER/${radius}/g" parameters.yml
    sed -i "s/RADIUS_PLACEHOLDER/${radius}/g" submit.slurm

    sbatch submit.slurm
    cd ..
done

    

