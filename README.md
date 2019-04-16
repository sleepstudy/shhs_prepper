# SHHS Prepper
ETL tool for the National Sleep Research Resource SHHS dataset

## How to use
To get configurable parameters available for you when using the script, simply run `python shhs-prepper.py --help`
```shell
usage: shhs-prepper.py [-h] [--batch_size BATCH_SIZE]
                       [--population_size POPULATION_SIZE] [--interspersed]
                       [--patient_ids PATIENT_IDS [PATIENT_IDS ...]]
                       [--path PATH] [--dump_path DUMP_PATH] [--profile]

Turns an shhs dataset into a useful dataframe/CSV set

optional arguments:
  -h, --help            show this help message and exit
  --batch_size BATCH_SIZE
                        How many EDF samples to batch together. Affects column
                        size
  --population_size POPULATION_SIZE
                        How many patients should be used?
  --interspersed        Should different EEG sensors be interleaved into the
                        same/different row?
  --path PATH           Path to SHHS folder
  --dump_path DUMP_PATH
                        Path to CSV result dump
  --profile             Profile mode?
```

## What does it do?
Upon running the script, it creates a CSV file and dumps it to a path (by default, set to `results.csv` on the current working directory. 
  * `batch_size` controls how many EDF samplings are recorded together as one CSV row. This means that given an EEG signal sampled
  at 125Hz (125 signal samplings per second), a `batch_size` of 125 would store 1 second of EEG data samples per row, with 125 columns 
  representing each sampled signal for a given EEG sensor (there are 2 on the SHHS dataset)
  * `interspersed` controls how many different EEG sensors should be interleaves into the same row. If passed (by running the program with 
  the `--interspersed` option), the 2 EEG sensors used in the SHHS dataset would be represented in the same row. If not, EEG sensors would 
  be represented into their own individual rows
  * `path` controls the path to the SHHS dataset folder. By default, it assumes the current working directory is the SHHS folder
  * `dump_path` controls the path and name of the resulting CSV file
  
## How does the output CSV look like?
  * With `batch_size` set to 2 and `--interspersed`, you would have the following columns per each row on your CSV:
    * `subject_id` - Patient Identifier
    * `in_cohort1` - Did the patient have a Cohort 1 EDF
    * `in_cohort2` - Did the patient have a Cohort 2 EDF
    * `cohort` - What cohort (1 or 2) was this EEG sampling drawn from
    * `sleep_stage` - What was the patient's sleep stage at this time?
    * `eeg1_0`, `eeg1_1`, `eeg2_0`, `eeg2_1` - With `--interleaved`, we store both `eeg1` and `eeg2` sensor data in the same row. 
    With `batch_size` set to 2, 2 EDF samplings for `eeg1` and `eeg2` are stored per row.
  * With `batch_size` set to 50, you would have the following columns per each row on your CSV:
    * `subject_id` - Patient Identifier
    * `in_cohort1` - Did the patient have a Cohort 1 EDF
    * `in_cohort2` - Did the patient have a Cohort 2 EDF
    * `cohort` - What cohort (1 or 2) was this EEG sampling drawn from
    * `sleep_stage` - What was the patient's sleep stage at this time?
    * `eeg_0`, `eeg_1` ... `eeg_49` - With `batch_size` set to 50, 50 EDF EEG samplings are stored per row.
    * `eeg_signal` - Tells us whether this represents the 1st or 2nd EEG signal (`1` or `2`), considering the 
    `--interspersed` option was not used.

## What's `--profile`
You most likely won't need to use it, but it runs `cProfiler` and to be honest, was an added option as I worked to improve the 
profiling of this script to run through across the SHHS monstrosity of a dataset.
