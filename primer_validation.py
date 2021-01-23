import csv
import Levenshtein
import pandas as pd
from IPython.display import display

primer_list = None
all_gene_list = None

output = None
analysis_seq_len =30


def primer_validation(_primer_path, _genelist_path, _output_path, _option='432', type='forward'):
    primer_list = get_primer_list(_primer_path)
    print(primer_list)
    all_gene_list = get_gene_list(_genelist_path, type)

    output = pd.DataFrame(data=[], index=range(0, len(all_gene_list)),
                          columns=create_initial_output_line(len(primer_list)))
    initialize_output(_genelist_path, output)

    for _row, _gene in enumerate(all_gene_list):
        _glen = len(_gene)

        for _col, _primer in enumerate(primer_list):
            _plen = len(_primer)
            _tmp_seq_distance = []
            for i in range(_glen - _plen + 1):
                _tmp = _gene[i:i + _plen]
                if _option == '321':
                    _tmp_seq_distance.append(validation_primer_321(_tmp, _primer))
                else:
                    _tmp_seq_distance.append(validation_primer_432(_tmp, _primer))

            if max(_tmp_seq_distance) == 0:
                output.iat[_row, _col + 2] = None
            else:
                output.iat[_row, _col + 2] = max(_tmp_seq_distance)
    print(_output_path)
    print(primer_list)
    # display(output)
    output.to_csv(r'' + _output_path, index=False, header=True)


def validation_primer_321(_tmp1, _tmp2):
    if (Levenshtein.distance(_tmp1, _tmp2) <= 1) & (_tmp1[-3:] == _tmp2[-3:]):
        return 1
    elif Levenshtein.distance(_tmp1, _tmp2) == 2:
        if (_tmp1[0] != _tmp2[0]) & (_tmp1[-3:] == _tmp2[-3:]):
            return 1
        else:
            return 0
    elif Levenshtein.distance(_tmp1, _tmp2) == 3:
        if (Levenshtein.distance(_tmp1[0:2], _tmp2[0:2]) == 3) & (_tmp1[-3:] == _tmp2[-3:]):
            return 1
        else:
            return 0
    else:
        return 0


def validation_primer_432(_tmp1, _tmp2):
    if (Levenshtein.distance(_tmp1, _tmp2) <= 2) & (_tmp1[-3:] == _tmp2[-3:]):
        return 1
    elif Levenshtein.distance(_tmp1, _tmp2) == 3:
        if (_tmp1[0] != _tmp2[0]) & (_tmp1[-3:] == _tmp2[-3:]):
            return 1
        else:
            return 0
    elif Levenshtein.distance(_tmp1, _tmp2) == 4:
        if (Levenshtein.distance(_tmp1[0:3], _tmp2[0:3]) == 4) & (_tmp1[-3:] == _tmp2[-3:]):
            return 1
        else:
            return 0
    else:
        return 0


def initialize_output(_path, _output):
    with open(_path, 'r') as file:
        reader = csv.reader(file)
        _seq_list = []
        _gene_list = []
        for row in reader:
            _seq_list.append((row[2][0:analysis_seq_len]))
            _gene_list.append((row[1][0:analysis_seq_len]))
        _seq_list = _seq_list[1:]
        _gene_list = _gene_list[1:]

    for _row in range(len(_seq_list)):
        _output.iat[_row, 0] = _seq_list[_row]
        _output.iat[_row, 1] = _gene_list[_row]


def create_initial_output_line(_iter):
    _tmp = ['seq', 'gene']
    for i in range(_iter):
        _tmp.append('primer_' + str(i))
    return _tmp


def get_primer_list(_path):
    with open(_path, 'r') as file:
        reader = csv.reader(file)
        _list = []
        for row in reader:
            _list.append((row[1]))
        return _list[1:]


def get_gene_list(_path, _type='forward'):
    with open(_path, 'r') as file:
        reader = csv.reader(file)
        _gene_list = []
        if _type == 'forward':
            for row in reader:
                _gene_list.append((row[2][0:analysis_seq_len]))
            return _gene_list[1:]
        else:
            for row in reader:
                row[2] = row[2][::-1]
                row[2] = row[2][1:]
                row[2] = make_complementary_seq(row[2])
                _gene_list.append((row[2][0:analysis_seq_len]))
            return _gene_list[1:]


def make_complementary_seq(_input_seq):
    output_seq = ''
    for tmp in _input_seq:
        if tmp == 'A':
            output_seq = output_seq + 'T'
        if tmp == 'G':
            output_seq = output_seq + 'C'
        if tmp == 'T':
            output_seq = output_seq + 'A'
        if tmp == 'C':
            output_seq = output_seq + 'G'
    return output_seq