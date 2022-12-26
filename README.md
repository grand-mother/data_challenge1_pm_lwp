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

## How to define an external pipeline

A pipeline is a chain of functions that pass results and parameters from one to another. In this case, we start with Efield traces and through a chain of operations want to arrive at ADC traces. That is what happens for one event. However, we can have multiple events in a file, and multiple files. Therefore, the execution of program is the following:

  1. Initialise the program
  2. Loop through the files
     a. For each file, read trees, initialise stuff needed at a file level, etc.
     b. Initialise an output file/files
  3. Loop through events inside a file
     a. Read the traces of the event
     b. Convert the traces to Voltage and ADC
     c. Fill the trees with traces
  4. After the event loop, write the file(s)
  
There are several operation blocks of which the pipeline may consist. They can belong to the different levels of the program described above:
1. Pre-file loop call -- a function called before the files loop. Only one per pipeline, **necessary**. Defines initial variables, such as necessary files paths, etc.
     a. name: the name of the function to be called
     b. type: prefileloop_call
     c. kwargs: parameter_name: parameter_value pairs that are passed to the function
     d. module: Python module where the function is defined
     e. return: a dictionary with variable_name: variable_value to be passed to the following elements of the pipeline. Can be empty.
2. Pre-event loop call -- a function called before the events loop. Useful for initialising data common for all events in the file, like filters coefficient, etc., that are not common for all files (due to, for example, possible different time bin for each run)
     a. name: the name of the function to be called
     b. type: preeventloop_call
     c. the rest the same as in the prefile-loop call
3. Call -- a function called inside the events loop.
     a. name: the name of the function to be called
     b. type: call
     c. the rest the same as in the prefile-loop call
4. add -- adds (with addition ;) ) an array to the array of the traces passed between the block
     a. name: addend - the name of the array to be added (must already be in the dictionary passed between blocks)
     b. type: add
     c. remarks: if a pair of [time_domain, frequency_domain] arrays is given, both are used. If only time_domain array is given, the traces in frequency domain are computed after the addition with (r)fft and returned with time traces.
5. add_randomized -- like add, but shuffles the addend array on the first index, that is assumed to be the traces (antennas) numbers
6. multiply -- multiplies the array of traces by the specified array
     a. name: multipland -- the name of the array by which the traces will be multiplied (must already be in the dictionary passed between blocks)
     b. type: multiply
     c. remarks: expects an array in frequency domain. Multiplies traces in frequency domain, then computes in time domain with i(r)fft
7. store -- stores the traces at this stage of pipeline to a tree
     a. name: not used, so can be anything
     b. type: store
     c. tree_type: the type of the tree in which to store the traces, such as VoltageEventTree, ADCEventTree, or in some strange cases, EfieldEventTree
     d. filename_suffix: the suffix that is added to the file name with Efield traces to create the file name of the file in which these traces will be stored. For example, for suffix "_voltage" we get Coarse2.root -> Coarse2_voltage.root.
  
Upon request, other blocks could be added, for example draw, for drawing traces at a specific point of the pipeline, or print, for printing traces.

All the data exchanged in the pipeline is stored in one dictionary. Traces in both time and frequency domains, as well as all the necessary coefficients (arrays), etc., are given to each block inside this dictionary. Each block should return a dictionary with variables (coefficient, traces, etc.) that it wants to pass on. The pipeline takes care abour updating the general dictionary with returned dictionary, and passing the general dictionary to each block.

The external pipeline files are stored in YAML format. For an example on how they should like, please look at PM_pipeline.yaml and XDU_pipeline.yaml. XDU pipeline does (rather unnecessarily) some stuff for each event, which PM pipeline does before each event loop - namely reading of the coefficients.

## How to define an internal pipeline

The logic is exactly the same as for the external pipeline. However, in this case, the whole pipeline is defined as a Python dictionary and this dictionary should be passed in your main function to "execute_pipeline()" function. An external YAML file can be translated to a Python dictionary, and that is exactly what happens when calling an external pipeline file (`yaml.safe_load("filename.yaml")`). Also, the python dictionary can be translated to a YAML file with a pyyaml call ( `yaml.safe_dump(pipeline_dict, sort_keys=False, stream=open("filename.yaml", "w")`).
