import gget 
import re
from Bio import Entrez 
import json
import random
import pandas as pd
import numpy as np

# input 
alpl_sequence = 'ATGATTTCACCATTCTTAGTACTGGCCATTGGCACCTGCCTTACTAACTCCTTAGTGCCAGAGAAAGAGAAAGACCCCAAGTACTGGCGAGACCAAGCGCAAGAGACACTGAAATATGCCCTGGAGCTTCAGAAGCTCAACACCAACGTGGCTAAGAATGTCATCATGTTCCTGGGAGATGGGATGGGTGTCTCCACAGTGACGGCTGCCCGCATCCTCAAGGGTCAGCTCCACCACAACCCTGGGGAGGAGACCAGGCTGGAGATGGACAAGTTCCCCTTCGTGGCCCTCTCCAAGACGTACAACACCAATGCCCAGGTCCCTGACAGTGCCGGCACCGCCACCGCCTACCTGTGTGGGGTGAAGGCCAATGAGGGCACCGTGGGGGTAAGCGCAGCCACTGAGCGTTCCCGGTGCAACACCACCCAGGGGAACGAGGTCACCTCCATCCTGCGCTGGGCCAAGGACGCTGGGAAATCTGTGGGCATTGTGACCACCACGAGAGTGAACCATGCCACCCCCAGCGCCGCCTACGCCCACTCGGCTGACCGGGACTGGTACTCAGACAACGAGATGCCCCCTGAGGCCTTGAGCCAGGGCTGTAAGGACATCGCCTACCAGCTCATGCATAACATCAGGGACATTGACGTGATCATGGGGGGTGGCCGGAAATACATGTACCCCAAGAATAAAACTGATGTGGAGTATGAGAGTGACGAGAAAGCCAGGGGCACGAGGCTGGACGGCCTGGACCTCGTTGACACCTGGAAGAGCTTCAAACCGAGATACAAGCACTCCCACTTCATCTGGAACCGCACGGAACTCCTGACCCTTGACCCCCACAATGTGGACTACCTATTGGGTCTCTTCGAGCCAGGGGACATGCAGTACGAGCTGAACAGGAACAACGTGACGGACCCGTCACTCTCCGAGATGGTGGTGGTGGCCATCCAGATCCTGCGGAAGAACCCCAAAGGCTTCTTCTTGCTGGTGGAAGGAGGCAGAATTGACCACGGGCACCATGAAGGAAAAGCCAAGCAGGCCCTGCATGAGGCGGTGGAGATGGACCGGGCCATCGGGCAGGCAGGCAGCTTGACCTCCTCGGAAGACACTCTGACCGTGGTCACTGCGGACCATTCCCACGTCTTCACATTTGGTGGATACACCCCCCGTGGCAACTCTATCTTTGGTCTGGCCCCCATGCTGAGTGACACAGACAAGAAGCCCTTCACTGCCATCCTGTATGGCAATGGGCCTGGCTACAAGGTGGTGGGCGGTGAACGAGAGAATGTCTCCATGGTGGACTATGCTCACAACAACTACCAGGCGCAGTCTGCTGTGCCCCTGCGCCACGAGACCCACGGCGGGGAGGACGTGGCCGTCTTCTCCAAGGGCCCCATGGCGCACCTGCTGCACGGCGTCCACGAGCAGAACTACGTCCCCCACGTGATGGCGTATGCAGCCTGCATCGGGGCCAACCTCGGCCACTGTGCTCCTGCCAGCTCGGCAGGCAGCCTTGCTGCAGGCCCCCTGCTGCTCGCGCTGGCCCTCTACCCCCTGAGCGTCCTGTTCTGA'

# HF overhangs
overhangs_df = pd.read_csv('HF_4bpOverhangsCombinations.csv', sep='\t')
overhang_setfiltered = overhangs_df[overhangs_df['Number of overhangs'] == 30]
overhang_set30 = overhang_setfiltered['Individual overhang sequences '].values.tolist()

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

sequence_overhang_dic = optimal_overhangs_search(alpl_sequence, overhang_set30)


def find_optimal_split_points(sequence: str, overhang_set: list, min_length: int, max_length: int):

    sequence_overhang_dic = optimal_overhangs_search(sequence, overhang_set)
    sequence_overhang_df = pd.DataFrame.from_dict(sequence_overhang_dic)

    max_number_of_parts = len(sequence)//min_length
    min_number_of_parts = len(sequence)//max_length
    print(f'max number of parts: {max_number_of_parts}')
    print(f'min number of parts: {min_number_of_parts}')

    # 3,4, or 5 parts --> cut sequence part - 1 
    part_dic_indexes = {}
    for split in range(min_number_of_parts,max_number_of_parts+1):
        print(f'n: {split}')
        
        end_indexes = [0]
        for part in range(1,split):
            print(f'part: {part}')
            
            # Randomly choose two indexes from sequence overhang dic to use as end delimeters
            valid_end_indexes = sequence_overhang_df[(sequence_overhang_df['end_index'] < (part*max_length)) & (sequence_overhang_df['end_index'] > (part*min_length))]
            print(f'valid_end_indexes: {valid_end_indexes}')

            random_end_index = random.choice(valid_end_indexes.end_index.values.tolist())
            print(f'random_end_index: {random_end_index}')

            end_indexes.append(random_end_index)

        end_indexes.append(len(sequence))

        part_dic_indexes[split] = {
            'end_indexes': sorted(end_indexes)
        }

    print(part_dic_indexes)
    return part_dic_indexes


def split_sequence(sequence: str, split_indexes: list):
    pass

