# bbi-sge
    The bioinformatic pipeline for Saturation Genome Editing in Brotman Baty Insitute

## News
    0.9.5 - update with csv read

## Documentation
    https://bbi.atlassian.net/wiki/spaces/SGE/pages/274202632/SGE+Pipeline

## Installation

1. You need to git clone follwing:  

        module load git/2.18.0
        git clone https://github.com/bbi-lab/bbi-sge.git

2. Download nextflow inside of bbi-sge folder that is created by step 1.   
        curl -fsSL https://get.nextflow.io | bash

3. Copy Seqprep folder to inside of bbi-sge folder
    scp -r uw-id@/net/gs/vol1/home/gf2/bin/SeqPrep/ uw-id@the directory of your nextflow

## Useage

How to run the pipeline:
1. Logining grid-head
2. Getting in tmux
3. Getting in qlogin
4. Changing parameters that you want in params.config
5. run the nextflow

Detail description of each steps

1. Login to gs cluster  

2. As the Nextflow pipeline is run interactively, please use a terminal multiplexer such as tmux. tmux sessions are persistent which means that programs in tmux will continue to run even if you get disconnected. You can start a tmux session by using:
    module load libevent/2.1.8 tmux/2.8
    tmux   
If you get disconnected, you can return to the head node you were using (grid-head1 or grid-head2) and type:
    tmux attach

3. Always start with a qlogin session before you begin the pipeline. This can be done using:
    qlogin -l mfree=48G 
You can submit the job without qlogin, but if the cluster is running out of rescource(memory), it will crush a middle of run.

4. Now, you can look up the params.config file to change parameters (you can change name of params.config file whatever you want it). Information about what each parameters for are descripted in config file.

5. Running the nextflow
        ./nextflow run main.nf -c <config name>.config 

All the required modules are loaded within the nextflow.config and main.nf, so you do not need to worry about loading modules 

General gs-cluster modules 
Required to load other modules; .bashrc
    module load modules modules-init modules-gs  

## Extra Features
Those are implanted features on nextflow that can be useful to run the pipeline.

1. resume  
This will resume the pipeline from where the problem arise after fixing the problem.
    nextflow run bbi-sge -c <config name> -resume

2. notification  
You can put tagging on the command line to get a notification when you are done
    nextflow run bbi-sge -c <config name> -N <recipient address>  

## dependency
    The script rquires Nextflow version >=21.10.6  

    gmp/6.1.2 mpfr/4.0.1 mpc/1.1.0 gcc/8.1.0 bcl2fastq/2.20 java/1.8.0 fastqc/0.11.7 EMBOSS/6.6.0
