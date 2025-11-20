# Project 3 Nextflow Template

For this project, remember to keep in a few things:

1. Most of the required references and files can be found in your `nextflow.config`

2. Make sure you give each process a label to request an appropriate amount of resources

3. Use the singularity containers provided on the website directions for the project

4. I have given you valid stub commands that will let you troubleshoot your workflow logic using the `-stub-run` command
- The stub-run commands assume that the first element in the tuple from the initial channel is named `sample_id` in processes
- Ensure that the appropriate inputs for certain processes are a tuple with the first element being the name from the initial channel
- The findPeaks stub will not be the same as `sample_id`. Remember that you will need to run findPeaks using the paired samples
(IP_rep1 + INPUT_rep1) and (IP_rep2 + INPUT_rep2). You should name the peak outputs using the replicate (i.e. rep1_peaks.txt and rep2_peaks.txt)
- You may alter the names used in the stub-run if it's easier for you

The stub runs assume that you have something like below so that it can name the fake files using the sample names - this will ensure
that your stub runs execute the same number of processes as the full pipeline should.
```
input:
tuple val(sample_id), path(file)
```

5. Use the subsampled data to start out with - you may need to eventually switch to the full data before your
pipeline is technically complete as sometimes peak calling may fail if not given enough input reads. 
- When the pipeline is working, change the `params` value in the original channel to the params encoding the
location of the full_samplesheet.csv

6. To remove regions using the blacklist, there are optional flags available in the `bedtools intersect` command

7. Create a single jupyter notebook that contains all of the results / figures and your write-up