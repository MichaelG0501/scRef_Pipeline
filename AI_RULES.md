**HPC EXECUTION & FILE SAFETY PROTOCOLS**

**0. Core Paths & Envs**
* **WD (Save Dir)**: `/rds/general/project/tumourheterogeneity1/ephemeral/VS_Code` (All created files MUST be saved here).
* **Data Dir**: Will be specified by the user per task.
* **Conda Init**: `eval "$(~/miniforge3/bin/conda shell.bash hook)"`
* **Env R (Default)**: `source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp`
* **Env R (UCell/GeneNMF)**: `source activate /rds/general/user/sg3723/home/anaconda3/envs/gnmf`
* **Env Python (Default)**: `source activate bidcell_temp`
* *Note*: Create new conda env if new software requires it.

**1. Execution Routing**
* **Interactive First (<= 8 cores / 64GB)**: For light/visualization tasks, generate ONLY the `.R` or `.py` script. Do NOT create `.sh` wrappers. The user will execute these manually in an interactive IDE (e.g., RStudio).
* **PBS Jobs (qsub)**: Mandatory for tasks exceeding limits, heavy parallelism, intensive I/O, or anything risky.

**2. PBS Job Script Standards**
* **Resource Sizing**: Dynamically determine optimal `ncpus`, `mem`, and `walltime`. Pick lowest safe values to minimize queue time while preventing OOM/timeout kills.
* **Single Job Template**:
    ```bash
    #!/bin/bash
    #PBS -l select=1:ncpus=<OPT_CPUS>:mem=<OPT_MEM>gb
    #PBS -l walltime=<OPT_TIME>
    #PBS -N <jobname>
    echo $(date +%T)
    module purge
    module load tools/dev
    <Conda Init>
    <Env>
    WD=/rds/general/project/tumourheterogeneity1/ephemeral/Auto_AG
    cd $WD
    <command>
    echo $(date +%T)
    ```
* **Parallel Jobs (`_master.sh`)**: If parallelism helps, write a master script to iterate folders and throttle submissions to max 46:
    ```bash
    for sample_folder in <SPECIFIED_DATA_DIR>/*_*_*/; do
      while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]; do sleep 180; done
      sample=$(basename "$sample_folder")
      qsub -v sample=$sample -N <real_script>.sh
    done
    ```

**3. File Modification & Permissions**
* **NO Permission Needed**:
    * **Creating files**: Must prefix persistent files with `Auto_` (e.g., `Auto_run.R`).
    * **Adding lines to existing files**: New code MUST be entirely wrapped within 20-hash comment blocks:
        `####################`
        `[your new code here]`
        `####################`
* **ASK Permission**:
    * Deleting/removing existing files.
    * Modifying or deleting existing lines outside the 20-hash blocks.

**4. Test Scripts & Cleanup**
* **NO Permission Needed** to create/run quick check scripts (data structure, packages).
* Always include `<Conda Init>` and appropriate `<Env>` before running.
* **Mandatory**: Delete test file immediately after run. 
* If uncertain it will delete, name it `delete_<desc>.R/.py`.
