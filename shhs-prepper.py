import argparse
import os.path
import itertools
from collections import namedtuple
from xml.etree import ElementTree
import math
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
  
  def get_sleepstages(edf_file):
    """
    Given an EDF file, returns a list of sleep stage values, 
    recorded in 30s intervals
    """
    return [
      sleep_stage_elem.text
      for sleep_stage_elem in ElementTree.parse(
        edf_file
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
  # [(sensor1_val, sensor2_val)]
  eeg_cohort_1_data = []
  eeg_cohort_1_sampling_frequency = (0, 0)
  eeg_cohort_1_present = False
    
  # [(sensor1_val, sensor2_val)]
  eeg_cohort_2_data = []
  eeg_cohort_2_sampling_frequency = (0, 0)
  eeg_cohort_2_present = False 
    
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
    Returns ((sensor1, sensor2), (freq1, freq2), [(sensor1, sensor2)], (freq1, freq2))
    """
    cohort_1 = []
    cohort_1_freq = (0, 0)
    cohort_2 = []
    cohort_2_freq = (0, 0)
    
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
        cohort_1 = list(zip(
          reader_1.readSignal(2).tolist(), 
          reader_1.readSignal(7).tolist()
        ))
        cohort_1_freq = (
          reader_1.getSampleFrequency(2),
          reader_1.getSampleFrequency(7)
        )
      
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
        cohort_2 = list(zip(
          reader_2.readSignal(2).tolist(), 
          reader_2.readSignal(7).tolist()
        ))
        cohort_2_freq = (
          reader_2.getSampleFrequency(2),
          reader_2.getSampleFrequency(7)
        )
      
    return cohort_1, cohort_1_freq, cohort_2, cohort_2_freq
  
  def generator():
    for patient_info in patient_info_generator():
      if patient_info.subject_id not in eeg_subject_ids:
        continue

      patient_info.__class__= PatientEEGInfo
      
      cohort1_data, cohort1_sf, cohort2_data, cohort2_sf = eeg_from_edfreader(
        patient_info.subject_id
      )
      patient_info.eeg_cohort_1_data = cohort1_data
      patient_info.eeg_cohort_1_sampling_frequency = cohort1_sf
      patient_info.eeg_cohort_2_data = cohort2_data
      patient_info.eeg_cohort_2_sampling_frequency = cohort2_sf
      
      patient_info.eeg_cohort_1_present = patient_info.subject_id in subject_ids_cohort_1
      patient_info.eeg_cohort_2_present = patient_info.subject_id in subject_ids_cohort_2

      yield patient_info
  
  return generator, eeg_subject_ids


class Settings:
  batch_size = 125
  population_size = 200
  interspersed = True
  patient_ids = []
  
  
def transform_to_dict(patient_eeg_info_generator, settings):
  
  def eeg_generator(eeg_data, sleep_stages_data, interspersed):
    eeg_data_len = len(eeg_data)
    sleepstage_len = len(sleep_stages_data)
    
    sleepstage_per_sensor = sleepstage_len/ float(
      eeg_data_len
    ) if float(eeg_data_len) > 0 else 0
    
    for eeg_index in range(0, eeg_data_len, settings.batch_size):
      eeg1_data, eeg2_data = zip(*eeg_data)
      eeg1_data_tuples = [
        ("eeg1_{0}".format(i), eeg)
        for i, eeg, in enumerate(eeg1_data[eeg_index:eeg_index+settings.batch_size])
      ]
      eeg2_data_tuples = [
        ("eeg2_{0}".format(i), eeg)
        for i, eeg, in enumerate(eeg2_data[eeg_index:eeg_index+settings.batch_size])
      ]
      sleep_stage = [
        (
          "sleep_stage", 
          sleep_stages_data[
            int(math.floor(i * sleepstage_per_sensor))
          ]
        )
      ]
      
      if interspersed:
        yield eeg1_data_tuples + eeg2_data_tuples + sleep_stage
      else:
        yield eeg1_data_tuples + sleep_stage
        yield eeg2_data_tuples + sleep_stage
  
  def eeg_info_transformer(patient_eeg_info):
    iterables = []
    dict_addons = [
      ("subject_id", patient_eeg_info.subject_id),
      ("in_cohort1", patient_eeg_info.eeg_cohort_1_present),
      ("in_cohort2", patient_eeg_info.eeg_cohort_2_present)
    ]
    
    if patient_eeg_info.eeg_cohort_1_present:
      dict_addons.append(("cohort", "1"))
      
      iterables.append(itertools.imap(lambda x: dict(x + dict_addons), eeg_generator(
        patient_eeg_info.eeg_cohort_1_data, 
        patient_eeg_info.sleep_stages_cohort_1, 
        settings.interspersed
      )))
      
    elif patient_eeg_info.eeg_cohort_2_present:
      dict_addons.append(("cohort", "2"))
      
      iterables.append(itertools.imap(lambda x: dict(x + dict_addons), eeg_generator(
        patient_eeg_info.eeg_cohort_2_data, 
        patient_eeg_info.sleep_stages_cohort_2, 
        settings.interspersed
      )))
      
    return itertools.chain.from_iterable(iterables)
  
  return itertools.chain.from_iterable(
    itertools.imap(eeg_info_transformer, patient_eeg_info_generator())
  )


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
  
  parsed = parser.parse_args()
  
  settings = Settings()
  shhs_dir_path = os.path.abspath(parsed.path)
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
  
  
  patient_info_generator, subject_ids = itertools.islice(
    get_patient_eeg(
      shhs_dir_path, 
      *get_patient_info(shhs_dir_path)
    ), 
    parsed.population_size
  )
    
  for patient_info in transform_to_dict(patient_info_generator, settings):
    print(len(patient_info.keys()))

  
main()