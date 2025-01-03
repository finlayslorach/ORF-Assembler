import gget 
import re
from Bio import Entrez 
import json
import random
import pandas as pd
import numpy as np

# input 
alpl_sequence = 'ATGATTTCACCATTCTTAGTACTGGCCATTGGCACCTGCCTTACTAACTCCTTAGTGCCAGAGAAAGAGAAAGACCCCAAGTACTGGCGAGACCAAGCGCAAGAGACACTGAAATATGCCCTGGAGCTTCAGAAGCTCAACACCAACGTGGCTAAGAATGTCATCATGTTCCTGGGAGATGGGATGGGTGTCTCCACAGTGACGGCTGCCCGCATCCTCAAGGGTCAGCTCCACCACAACCCTGGGGAGGAGACCAGGCTGGAGATGGACAAGTTCCCCTTCGTGGCCCTCTCCAAGACGTACAACACCAATGCCCAGGTCCCTGACAGTGCCGGCACCGCCACCGCCTACCTGTGTGGGGTGAAGGCCAATGAGGGCACCGTGGGGGTAAGCGCAGCCACTGAGCGTTCCCGGTGCAACACCACCCAGGGGAACGAGGTCACCTCCATCCTGCGCTGGGCCAAGGACGCTGGGAAATCTGTGGGCATTGTGACCACCACGAGAGTGAACCATGCCACCCCCAGCGCCGCCTACGCCCACTCGGCTGACCGGGACTGGTACTCAGACAACGAGATGCCCCCTGAGGCCTTGAGCCAGGGCTGTAAGGACATCGCCTACCAGCTCATGCATAACATCAGGGACATTGACGTGATCATGGGGGGTGGCCGGAAATACATGTACCCCAAGAATAAAACTGATGTGGAGTATGAGAGTGACGAGAAAGCCAGGGGCACGAGGCTGGACGGCCTGGACCTCGTTGACACCTGGAAGAGCTTCAAACCGAGATACAAGCACTCCCACTTCATCTGGAACCGCACGGAACTCCTGACCCTTGACCCCCACAATGTGGACTACCTATTGGGTCTCTTCGAGCCAGGGGACATGCAGTACGAGCTGAACAGGAACAACGTGACGGACCCGTCACTCTCCGAGATGGTGGTGGTGGCCATCCAGATCCTGCGGAAGAACCCCAAAGGCTTCTTCTTGCTGGTGGAAGGAGGCAGAATTGACCACGGGCACCATGAAGGAAAAGCCAAGCAGGCCCTGCATGAGGCGGTGGAGATGGACCGGGCCATCGGGCAGGCAGGCAGCTTGACCTCCTCGGAAGACACTCTGACCGTGGTCACTGCGGACCATTCCCACGTCTTCACATTTGGTGGATACACCCCCCGTGGCAACTCTATCTTTGGTCTGGCCCCCATGCTGAGTGACACAGACAAGAAGCCCTTCACTGCCATCCTGTATGGCAATGGGCCTGGCTACAAGGTGGTGGGCGGTGAACGAGAGAATGTCTCCATGGTGGACTATGCTCACAACAACTACCAGGCGCAGTCTGCTGTGCCCCTGCGCCACGAGACCCACGGCGGGGAGGACGTGGCCGTCTTCTCCAAGGGCCCCATGGCGCACCTGCTGCACGGCGTCCACGAGCAGAACTACGTCCCCCACGTGATGGCGTATGCAGCCTGCATCGGGGCCAACCTCGGCCACTGTGCTCCTGCCAGCTCGGCAGGCAGCCTTGCTGCAGGCCCCCTGCTGCTCGCGCTGGCCCTCTACCCCCTGAGCGTCCTGTTCTGA'

with open('config.json', 'r') as f: config = json.load(f)
overhang_set_30_file = config['OVERHANG_SET_30']
with open(overhang_set_30_file, 'r') as f: overhang_set_30 = f.read()
print("overhang_set_30:")
print(overhang_set_30)

# Are optimal overhangs in the sequence 
def optimal_overhangs_search(sequence: str, overhang_set: list):
    overhang_set_str = re.sub(r'[\[\]\s\']+', '', str(overhang_set))
    pattern = f"({overhang_set_str.replace(',', '|')})"
    sequence_overhang_dic = {'overhang': [], 'start_index': [], 'end_index': []}
    for match in re.finditer(pattern, str(sequence.strip())):
        sequence_overhang_dic['overhang'].append(match.group())
        sequence_overhang_dic['start_index'].append(match.start())
        sequence_overhang_dic['end_index'].append(match.end())
    return sequence_overhang_dic


def find_optimal_split_points(sequence: str, min_number_of_parts: int, max_number_of_parts: int, max_length: int,
                            min_length: int, sequence_overhang_dic: dict):
    sequence_overhang_df = pd.DataFrame.from_dict(sequence_overhang_dic)
    optimal_split_points = {}
    for split in range(min_number_of_parts,max_number_of_parts+1):
        end_indexes = [0]
        for part in range(1,split):            
            # Randomly choose two indexes from sequence overhang dic to use as end delimeters
            valid_end_indexes = sequence_overhang_df[(sequence_overhang_df['end_index'] < (part*max_length)) & (sequence_overhang_df['end_index'] > (part*min_length))]
            random_end_index = random.choice(valid_end_indexes.end_index.values.tolist())
            end_indexes.append(random_end_index)
        end_indexes.append(len(sequence))
        optimal_split_points[split] = {'end_indexes': sorted(end_indexes)}
    return optimal_split_points
    

def main(sequence: str, min_length: int, max_length: int, overhang_set: list):
    ## User parameters for splitting sequence
    max_number_of_parts = len(sequence)//min_length
    min_number_of_parts = len(sequence)//max_length
    print(f'max number of parts: {max_number_of_parts}')
    print(f'min number of parts: {min_number_of_parts}')

    ## Find all optimal overhangs in the sequence
    sequence_overhang_dic = optimal_overhangs_search(sequence, overhang_set)
    sequence_overhang_df = pd.DataFrame.from_dict(sequence_overhang_dic)

    ## Find a set of optimal split points
    optimal_split_points = {}
    for split in range(min_number_of_parts,max_number_of_parts+1):        
        end_indexes = [0]
        for part in range(1,split):
            # Randomly choose two indexes from sequence overhang dic to use as end delimeters
            valid_end_indexes = sequence_overhang_df[(sequence_overhang_df['end_index'] < (part*max_length)) & (sequence_overhang_df['end_index'] > (part*min_length))]
            random_end_index = random.choice(valid_end_indexes.end_index.values.tolist())
            end_indexes.append(random_end_index)
        end_indexes.append(len(sequence))
        optimal_split_points[split] = {'end_indexes': sorted(end_indexes)}

    ## Split the sequence into parts
    split_sequence_options = {}
    for split_type, split_index_list in optimal_split_points.items():
        split_sequence_list = []
        for index in range(0, len(split_index_list['end_indexes'])-1):
            split_sequence = sequence[split_index_list['end_indexes'][index]:split_index_list['end_indexes'][index+1]]
            split_sequence_list.append(split_sequence) 
        split_sequence_options[split_type] = split_sequence_list
    return split_sequence_options

main(alpl_sequence, 300, 500, overhang_set_30)



     
    

