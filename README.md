# A contribution to the Data Challenge 1
*by Pragati Mitra and Lech Wiktor Piotrowski*

## Prerequisites (necessary!):
1. Please use the root_trees.py (`grand.io.root_trees`) from the latest `dev_io_root` branch of grandlib. If you can't get it, I can send you a zip and try to help you use it.
2. For PM pipeline, please download the necessary PM files (antenna model, coefficients, etc.). 300 MB :(
   1. temporarily at: https://mega.nz/file/YZZXxQAZ#imMdiT8ficDi2XIMaHSu56ZzcuD61S0Uauj3pkUkX6o
   2. unpack into `PM_functions` (so the structure is `data_challenge1_pm_lwp/PM_functions/PM_files`)

      (execute `tar -xvzf PM_files.tar.gz` in `data_challenge1_pm_lwp/PM_functions`)
3. For XDU pipeline, please download the necessary XDU files (antenna model, coefficients, etc.). 1 GB :(
   1. temporarily at: https://mega.nz/file/kFhDGTzK#UCuP_usSD4TxdxsmtbKSd_Z_1dP66MFFprrXIuJ6BdA
   2. unpack into `XDU_electronic_chain` (so the structure is `data_challenge1_pm_lwp/XDU_electronic_chain/XDU_files`)

        (execute `tar -xvzf XDU_files.tar.gz` in `data_challenge1_pm_lwp/XDU_electronic_chain`)
4. Please install yaml python module
   
   `pip install yaml` as root should do.


## Remarks:
1. I assume that the source root files are in the `data_challenge1_pm_lwp/data/` subdirectory (but you can change, of course)
2. I copied the data_challenge1 Coarse*.root files into `data_challenge1_pm_lwp/data/`
3. I store into `results/` subdirectory (but you can change)
4. Before each re-run the `results/` subdirectory contents need to be deleted manually (easiest on command line: `rm -rf results/; ./thenameofthescript.py someparameters`)
5. I'll add an option to overwrite the files later
6. If you want repeatable clean results, comment out the galactic noise addition

## How to run ###

The pipeline can be run in **several ways**:
1. With **external pipeline config file** *(most universal)*:
  
  * Executes the pipeline specified in PM_pipeline.yaml file

     `rm -rf results; ./efield2adc_alltraces_external_pipeline.py data/*.root -p PM_pipeline.yaml -od results`

    (replace PM with XDU for XDU pipeline)

2. With **internal pipeline config file**:

  * Easy:
    
    Executes the pipeline internally stored in the called script as a dictionary
 
     `rm -rf results; ./efield2adc_alltraces_internal_pipeline_XDU.py data/*.root -od results/`

    (replace XDU with PM for PM pipeline)

  * For full control:
  
    The script calls the pipeline specified programmatically (python code) in `XDU_electronic_chain/XDU_manual_pipeline.py`

       `rm -rf results; ./efield2adc_alltraces_XDU.py data/*.root -od results`

    (replace XDU with PM for PM pipeline)

  * If trace by trace analysis is needed (much slower), also calls the pipeline specified in `XDU_electronic_chain/XDU_manual_pipeline.py`:

       `rm -rf results; ./efield2adc_tracebytrace_XDU.py data/*.root -od results`

    (replace XDU with PM for PM pipeline)  

As you probably guessed:
  1. First parameter is input file(s)
  2. `-od` selects output directory
  3. `-p` is the external pipline file for `efield2adc_alltraces_external_pipeline.py`
  
All command line parameters can be obtained with `--help` switch.

Currently, the two first options store voltage and ADC counts in separate files, while two last in same files. Can be modified.