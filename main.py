import os.path
from os import path
import csv
import Levenshtein
import numpy as np
import pandas as pd
from IPython.display import display
import primer_validation as pv

# input_path = './input/'
# init_human_primer_path = 'primersets/hphf_primer_list_final.csv'
# init_human_variable_path = 'human_heavy_variable.csv'
# output_path = './output/output_3213.csv'

input_path = './input/gene_list/mouse/'
input_gene_list_path = './input/gene_list/mouse/'
input_primer_candidates_path = './input/primer_candidates/mouse/'
output_path = './output_mouse/'


def main():
    # pv.primer_validation(input_path + init_human_primer_path, input_path + init_human_variable_path, output_path)

    # make_direc('./output/aa')

    make_output_direc(input_gene_list_path, output_path)

    direc = os.listdir(input_gene_list_path)
    for sub_direc in direc:
        if os.path.isdir(input_path + sub_direc):
            tmp_direc = os.listdir(input_path + sub_direc)
            for tmp in tmp_direc:
                if os.path.isdir(input_primer_candidates_path + sub_direc + '/' + tmp) == False:
                    continue
                if check_forward(tmp) == 'forward':
                    _primer_path = get_primer_csv_from_path(input_primer_candidates_path + sub_direc + '/' + tmp + '/')
                    _gene_path = get_genelist_csv_from_path(input_gene_list_path + sub_direc + '/' + tmp + '/')
                    _output_path = output_path + sub_direc + '/' + tmp + '/output_432.csv'
                    pv.primer_validation(_primer_path, _gene_path, _output_path)
                    _output_path = output_path + sub_direc + '/' + tmp + '/output_321.csv'
                    pv.primer_validation(_primer_path, _gene_path, _output_path, _option='321')
                if check_forward(tmp) == 'reverse':
                    _primer_path = get_primer_csv_from_path(input_primer_candidates_path + sub_direc + '/' + tmp + '/')
                    _gene_path = get_genelist_csv_from_path(input_gene_list_path + sub_direc + '/' + tmp + '/',
                                                            'reverse')
                    _output_path = output_path + sub_direc + '/' + tmp + '/output_432.csv'
                    pv.primer_validation(_primer_path, _gene_path, _output_path, type='reverse')
                    _output_path = output_path + sub_direc + '/' + tmp + '/output_321.csv'
                    pv.primer_validation(_primer_path, _gene_path, _output_path, _option='321', type='reverse')


def get_primer_csv_from_path(_path):
    direc = os.listdir(_path)
    for tmp in direc:
        if 'primer' in tmp:
            return _path + tmp
    return os.error()


def get_genelist_csv_from_path(_path, option='forward'):
    direc = os.listdir(_path)
    for tmp in direc:
        if option == 'forward':
            if 'variable' in tmp:
                return _path + tmp
        else:
            if 'joining' in tmp:
                return _path + tmp
    return _path


def check_forward(_path):
    if 'forward' in _path:
        return 'forward'
    else:
        return 'reverse'


def make_output_direc(_input_path, _output_path):
    make_direc(_output_path)
    direc = os.listdir(input_gene_list_path)
    for sub_direc in direc:
        if os.path.isdir(_input_path + sub_direc):
            make_direc(_output_path + sub_direc)
            tmp_direc = os.listdir(_input_path + sub_direc)
            print(tmp_direc)
            for tmp in tmp_direc:
                if os.path.isdir(_input_path + sub_direc + '/' + tmp):
                    make_direc(_output_path + sub_direc + '/' + tmp)


def make_direc(_path):
    if os.path.isdir(_path):
        print("directory %s is already exist" % _path)
    else:
        os.mkdir(_path)
        print("Successfully created the directory %s " % _path)


if __name__ == "__main__":
    # execute only if run as a script
    main()
