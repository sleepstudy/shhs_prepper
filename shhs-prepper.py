import argparse
import os.path
import itertools
from collections import namedtuple
from xml.etree import ElementTree
import math
import csv
import pyedflib


def get_matching_filenames(cohort_dir, prefix, suffix):
  """
  Returns a list of tuples including the name of the annotated PSG 
  and the normalized patient ID
  """
  return zip(*[
    (dir_entry, dir_entry.replace(prefix, '', 1).replace(suffix, '', 1))
    for dir_entry in os.listdir(cohort_dir)
    if dir_entry.startswith(prefix) and dir_entry.endswith(suffix)
  ])
  

class PatientInfo:
  subject_id = None
  cohort1 = False
  cohort2 = False
  
  # 30s intervals of sleep stages
  sleep_stages_cohort_1 = []
  sleep_stages_cohort_2 = []

  
def get_patient_info(shhs_dir_path):
  """
  Returns a generator of PatientInfo obtained from the SHHS dataset as 
  well as a list of subject IDs whose PatientInfo data is being returned
  """
  
  profusion_dir = os.path.normpath(
    os.path.join(shhs_dir_path, './polysomnography/annotations-events-profusion')
  )
  shhs1_dir = os.path.normpath(
    os.path.join(profusion_dir, './shhs1')
  )
  shhs2_dir = os.path.normpath(
    os.path.join(profusion_dir, './shhs2')
  )
  
  cohort1_xml_filenames, subject_ids_cohort_1 = get_matching_filenames(
    shhs1_dir, 
    prefix='shhs1-', 
    suffix='-profusion.xml'
  )
  cohort2_xml_filenames, subject_ids_cohort_2 = get_matching_filenames(
    shhs2_dir, 
    prefix='shhs2-', 
    suffix='-profusion.xml'
  )
  all_subject_ids = set(subject_ids_cohort_1 + subject_ids_cohort_2)
  
  patient_triplets = [(
    subject_id, 
    subject_id in subject_ids_cohort_1,
    subject_id in subject_ids_cohort_2
  ) for subject_id in all_subject_ids]
  
  def get_sleepstages(xml_file):
    """
    Given an XML file, returns a list of sleep stage values, 
    recorded in 30s intervals
    """
    return [
      sleep_stage_elem.text
      for sleep_stage_elem in ElementTree.parse(
        xml_file
      ).getroot().find('SleepStages')
    ]
  
  def generator():
    for (subject_id, cohort1_exists, cohort2_exists) in patient_triplets:
      patient_info = PatientInfo()
      patient_info.subject_id = subject_id
      patient_info.cohort1 = cohort1_exists
      patient_info.cohort2 = cohort2_exists

      patient_info.sleep_stages_cohort_1 = get_sleepstages(
        os.path.normpath(
          os.path.join(
            shhs1_dir, 
            cohort1_xml_filenames[
              subject_ids_cohort_1.index(subject_id)
            ]
          )
        )
      ) if cohort1_exists else []
      patient_info.sleep_stages_cohort_2 = get_sleepstages(
        os.path.normpath(
          os.path.join(
            shhs2_dir, 
            cohort2_xml_filenames[
              subject_ids_cohort_2.index(subject_id)
            ]
          )
        )
      ) if cohort2_exists else []

      yield patient_info
  
  return generator, all_subject_ids

    
class PatientEEGInfo(PatientInfo):
  # Iterator[(sensor1_val, sensor2_val)]
  eeg_cohort_1_data = iter([])
  eeg_cohort_1_sampling_frequency = (0, 0)
  eeg_cohort_1_present = False
  eeg_cohort_1_len = 0
    
  # Iterator[(sensor1_val, sensor2_val)]
  eeg_cohort_2_data = iter([])
  eeg_cohort_2_sampling_frequency = (0, 0)
  eeg_cohort_2_present = False
  eeg_cohort_2_len = 0
    
  def __str__(self):
    return 'EEG Info: {subject_id}'.format(
      subject_id=self.subject_id
    )
    
    
def get_patient_eeg(shhs_dir_path, patient_info_generator, info_subject_ids):
  """
  Returns a generator of EEG objects obtained from the SHHS dataset as 
  well as a list of subject IDs whose EEG data is being returned
  """
  
  edfs_dir = os.path.normpath(
    os.path.join(shhs_dir_path, './polysomnography/edfs')
  )
  shhs1_dir = os.path.normpath(
    os.path.join(edfs_dir, './shhs1')
  )
  shhs2_dir = os.path.normpath(
    os.path.join(edfs_dir, './shhs2')
  )
  
  cohort1_edf_filenames, subject_ids_cohort_1 = get_matching_filenames(
    shhs1_dir, 
    prefix='shhs1-', 
    suffix='.edf'
  )
  cohort2_edf_filenames, subject_ids_cohort_2 = get_matching_filenames(
    shhs2_dir, 
    prefix='shhs2-', 
    suffix='.edf'
  )
  eeg_subject_ids = set(
    subject_ids_cohort_1 + subject_ids_cohort_2
  ).intersection(info_subject_ids)
  
  def eeg_from_edfreader(subject_id):
    """
    Returns (
      [(sensor1, sensor2)], (freq1, freq2), len
      [(sensor1, sensor2)], (freq1, freq2), len
    )
    """
    cohort_1 = []
    cohort_1_freq = (0, 0)
    cohort_1_len = 0
    cohort_2 = []
    cohort_2_freq = (0, 0)
    cohort_2_len = 0
    
    if subject_id in subject_ids_cohort_1:
      with pyedflib.EdfReader(
        os.path.normpath(
          os.path.join(
            shhs1_dir, 
            cohort1_edf_filenames[
              subject_ids_cohort_1.index(subject_id)
            ]
          )
        )
      ) as reader_1:
        signal_1 = reader_1.readSignal(2)
        signal_2 = reader_1.readSignal(7)
        cohort_1 = itertools.izip(signal_1, signal_2)
        cohort_1_freq = (
          reader_1.getSampleFrequency(2),
          reader_1.getSampleFrequency(7)
        )
        cohort_1_len = min(signal_1.size, signal_2.size)
      
    if subject_id in subject_ids_cohort_2:
      with pyedflib.EdfReader(
        os.path.normpath(
          os.path.join(
            shhs2_dir, 
            cohort2_edf_filenames[
              subject_ids_cohort_2.index(subject_id)
            ]
          )
        )
      ) as reader_2:
        signal_1 = reader_2.readSignal(2)
        signal_2 = reader_2.readSignal(7)
        cohort_2 = itertools.izip(signal_1, signal_2)
        cohort_2_freq = (
          reader_2.getSampleFrequency(2),
          reader_2.getSampleFrequency(7)
        )
        cohort_2_len = min(signal_1.size, signal_2.size)
      
    return cohort_1, cohort_1_freq, cohort_1_len, cohort_2, cohort_2_freq, cohort_2_len
  
  def generator():
    for patient_info in patient_info_generator():
      if patient_info.subject_id not in eeg_subject_ids:
        continue

      patient_info.__class__= PatientEEGInfo
      
      cohort1_data, cohort1_sf, cohort1_len, cohort2_data, cohort2_sf, cohort2_len = eeg_from_edfreader(
        patient_info.subject_id
      )
      patient_info.eeg_cohort_1_data = cohort1_data
      patient_info.eeg_cohort_1_sampling_frequency = cohort1_sf
      patient_info.eeg_cohort_1_len = cohort1_len
      patient_info.eeg_cohort_2_data = cohort2_data
      patient_info.eeg_cohort_2_sampling_frequency = cohort2_sf
      patient_info.eeg_cohort_2_len = cohort2_len
      
      patient_info.eeg_cohort_1_present = patient_info.subject_id in subject_ids_cohort_1
      patient_info.eeg_cohort_2_present = patient_info.subject_id in subject_ids_cohort_2

      yield patient_info
  
  return generator, eeg_subject_ids


class Settings:
  batch_size = 125
  population_size = 200
  interspersed = True
  patient_ids = []
  dump_path = None
  
  
def transform_to_dicts(settings, patient_eeg_info_generator):
  
  def eeg_generator(eeg_data, sleep_stages_data, interspersed, batch_size, eeg_data_len):
    
    sleepstage_per_sensor = len(sleep_stages_data) / float(
      eeg_data_len
    ) if float(eeg_data_len) > 0 else 0
    
    eeg1_data, eeg2_data = itertools.tee(eeg_data, 2)
    
    for eeg_index in range(0, eeg_data_len, batch_size):
      sleep_stage = sleep_stages_data[
        int(math.floor(eeg_index * sleepstage_per_sensor))
      ]
      
      if interspersed:
        row_dict = {
          "eeg1_{0}".format(i): eeg
          for i, eeg in enumerate(itertools.islice(eeg_data, batch_size))
        }
        row_dict.update({
          "eeg2_{0}".format(i): eeg
          for i, eeg, in enumerate(itertools.islice(eeg_data, batch_size))
        })
        row_dict["sleep_stage"] = sleep_stage
        
        yield row_dict
      else:
        row_dict = {
          "eeg_{0}".format(i): eeg
          for i, eeg, in enumerate(itertools.islice(eeg_data, batch_size))
        }
        row_dict["eeg_signal"] = "1"
        row_dict["sleep_stage"] = sleep_stage
        yield row_dict.items()
        
        row_dict = {
          "eeg_{0}".format(i): eeg
          for i, eeg, in enumerate(itertools.islice(eeg_data, batch_size))
        }
        row_dict["eeg_signal"] = "2"
        row_dict["sleep_stage"] = sleep_stage
        yield row_dict.items()
  
  def eeg_info_transformer(patient_eeg_info):
    iterables = []
    dict_addons = [
      ("subject_id", patient_eeg_info.subject_id),
      ("in_cohort1", patient_eeg_info.eeg_cohort_1_present),
      ("in_cohort2", patient_eeg_info.eeg_cohort_2_present)
    ]
    
    if patient_eeg_info.eeg_cohort_1_present:
      iterables.append(itertools.imap(
        lambda x: dict(x + dict_addons + [("cohort", "1")]), 
        eeg_generator(
          patient_eeg_info.eeg_cohort_1_data, 
          patient_eeg_info.sleep_stages_cohort_1,
          settings.interspersed,
          settings.batch_size,
          patient_eeg_info.eeg_cohort_1_len
        )
      ))
      
    if patient_eeg_info.eeg_cohort_2_present:
      iterables.append(itertools.imap(
        lambda x: dict(x + dict_addons + [("cohort", "2")]), 
        eeg_generator(
          patient_eeg_info.eeg_cohort_2_data, 
          patient_eeg_info.sleep_stages_cohort_2, 
          settings.interspersed,
          settings.batch_size,
          patient_eeg_info.eeg_cohort_2_len
        )
      ))
      
    return itertools.chain.from_iterable(iterables)
  
  return itertools.chain.from_iterable(
    itertools.imap(eeg_info_transformer, patient_eeg_info_generator())
  )


def profile(shhs_dir_path, settings):
  import cProfile
  import pstats
  
  pr = cProfile.Profile()
  
  pr.enable()
  patient_info_generator, subject_ids = itertools.islice(
    get_patient_eeg(
      shhs_dir_path, 
      *get_patient_info(shhs_dir_path)
    ), 
    100
  )
  iterator = itertools.islice(
    transform_to_dicts(settings, patient_info_generator),
    4
  )
  iterator.next()
  iterator.next()
  iterator.next()
  iterator.next()
  pr.disable()
  
  pstats.Stats(pr).sort_stats('cumulative').print_stats()


def main():
  """
  Columns are: 
    subject_id
    in_cohort1 
    in_cohort2
    cohort
    [eeg1_1 ... eeg1_batch_size] 
    [eeg2_1 ... eeg2_batch_size]
    sleep_stage
  """
  
  parser = argparse.ArgumentParser(
    description='Turns an shhs dataset into a useful dataframe/CSV set'
  )
  parser.add_argument(
    '--batch_size', 
    type=int, 
    default=125, 
    help='How many EDF samples to batch together. Affects column size'
  )
  parser.add_argument(
    '--population_size', 
    type=int, 
    default=200, 
    help='How many patients should be used?'
  )
  parser.add_argument(
    '--interspersed', 
    default=False, 
    action='store_true',
    help='Should different EEG sensors be interleaved into the same/different row?'
  )
  parser.add_argument(
    '--patient_ids', 
    nargs='+', 
    help='List of Patient IDs to target'
  )
  parser.add_argument(
    '--path', 
    type=str, 
    default='.', 
    help='Path to SHHS folder'
  )
  parser.add_argument(
    '--dump_path', 
    type=str, 
    default='./results.csv', 
    help='Path to CSV result dump'
  )
  parser.add_argument(
    '--profile', 
    default=False, 
    action='store_true',
    help='Profile mode?'
  )
  
  parsed = parser.parse_args()
  
  settings = Settings()
  shhs_dir_path = os.path.abspath(parsed.path)
  dump_path = os.path.abspath(parsed.dump_path)
  try:
    batch_size = int(parsed.batch_size)
    settings.batch_size = batch_size
    if batch_size < 0:
      raise ValueError
  except ValueError:
    print("Batch size invalid, using default of 125")
  try:
    population_size = int(parsed.population_size)
    settings.population_size = population_size
    if population_size < 0:
      raise ValueError
  except ValueError:
    print("Population size invalid, using default of 200")
  settings.interspersed = parsed.interspersed
  settings.dump_path = dump_path
  
  if parsed.profile:
    profile(shhs_dir_path, settings)
    return
  
  patient_info_generator, subject_ids = itertools.islice(
    get_patient_eeg(
      shhs_dir_path, 
      *get_patient_info(shhs_dir_path)
    ), 
    parsed.population_size
  )
  
  column_names = ["subject_id", "in_cohort1", "in_cohort2", "cohort", "sleep_stage"]
  if settings.interspersed:
    column_names = column_names + [
      "eeg1_{0}".format(i)
      for i in range(settings.batch_size)
    ] + [
      "eeg2_{0}".format(i)
      for i in range(settings.batch_size)
    ]
  else:
    column_names = column_names + [
      "eeg_{0}".format(i)
      for i in range(settings.batch_size)
    ] + ["eeg_signal"]
  
  with open(settings.dump_path, 'w+') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=column_names)
    writer.writeheader()
    for patient_info in transform_to_dicts(settings, patient_info_generator):
      print("Inserted row into CSV")
      writer.writerow(patient_info)

  
main()