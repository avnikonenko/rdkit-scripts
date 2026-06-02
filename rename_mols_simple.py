#!/usr/bin/env python3

__author__ = 'Pavel Polishchuk'

import argparse
import gzip
import re
import sys


FIELD_RE = re.compile(r'^\s*>\s*<([^>]*)>')


def read_names(names_fname):

    names = {}

    with open(names_fname) as f:
        for i, line in enumerate(f, 1):
            line = line.rstrip('\r\n')
            if not line:
                continue
            tmp = line.split('\t')
            if len(tmp) != 2:
                raise ValueError(f'Line {i} in {names_fname} should contain exactly two tab-separated columns')
            if tmp[0] in names:
                raise ValueError(f'Duplicated molecule name "{tmp[0]}" in {names_fname} on line {i}')
            names[tmp[0]] = tmp[1]

    return names


def remove_field(lines, field_name):

    res = []
    i = 0

    while i < len(lines):
        m = FIELD_RE.match(lines[i])
        if m and m.group(1) == field_name:
            i += 1
            while i < len(lines) and lines[i].strip():
                i += 1
            if i < len(lines):
                i += 1
        else:
            res.append(lines[i])
            i += 1
            if m:
                while i < len(lines):
                    res.append(lines[i])
                    line = lines[i]
                    i += 1
                    if not line.strip():
                        break

    return res


def rename_record(lines, names, old_name_field):

    old_name = lines[0].rstrip('\r\n')
    new_name = names.get(old_name, old_name)

    lines[0] = f'{new_name}\n'
    lines = remove_field(lines, old_name_field)
    lines.extend([f'>  <{old_name_field}>\n', f'{old_name}\n', '\n'])

    return lines, old_name in names


def open_input(input_fname):

    if input_fname is None or input_fname == '/dev/stdin':
        return sys.stdin, False

    input_fname_l = input_fname.lower()
    if input_fname_l.endswith('.sdf.gz'):
        return gzip.open(input_fname, 'rt'), True
    if input_fname_l.endswith('.sdf'):
        return open(input_fname), True
    raise ValueError('Input file format is not supported. It can be only SDF or SDF.GZ.')


def open_output(output_fname):

    if output_fname is None or output_fname == '/dev/stdout':
        return sys.stdout, False
    return open(output_fname, 'wt'), True


def rename_molecules(input_fname, output_fname, names_fname, old_name_field):

    names = read_names(names_fname)
    fin, close_fin = open_input(input_fname)
    fout, close_fout = open_output(output_fname)

    record = []
    absent = 0
    total = 0

    try:
        for line in fin:
            if line.strip() == '$$$$':
                if record:
                    record, found = rename_record(record, names, old_name_field)
                    total += 1
                    if not found:
                        absent += 1
                    fout.writelines(record)
                fout.write('$$$$\n')
                record = []
            else:
                record.append(line)

        if record:
            record, found = rename_record(record, names, old_name_field)
            total += 1
            if not found:
                absent += 1
            fout.writelines(record)
            fout.write('$$$$\n')
    finally:
        if close_fin:
            fin.close()
        if close_fout:
            fout.close()

    sys.stderr.write(f'{absent} of {total} molecules had titles absent from names file\n')


def main():

    parser = argparse.ArgumentParser(description='Rename molecule titles in an SDF file using a two-column '
                                                 'tab-separated text file with old and new names.')
    parser.add_argument('-i', '--input', metavar='FILENAME', required=False, default=None,
                        help='input SDF or SDF.GZ file. If omitted STDIN will be read as SDF.')
    parser.add_argument('-o', '--output', metavar='FILENAME', required=False, default=None,
                        help='output SDF file. If omitted output will be redirected to STDOUT.')
    parser.add_argument('-n', '--names', metavar='FILENAME', required=True,
                        help='tab-separated text file with old molecule titles in the first column and new titles in '
                             'the second column.')
    parser.add_argument('-f', '--old-name-field', metavar='STRING', required=False, default='_old_name',
                        help='field name where old molecule titles will be stored. Existing fields with this name '
                             'will be replaced silently. Default: _old_name.')
    args = parser.parse_args()

    if not args.old_name_field:
        parser.error('--old-name-field cannot be empty')
    if '<' in args.old_name_field or '>' in args.old_name_field:
        parser.error('--old-name-field cannot contain angle brackets')

    try:
        rename_molecules(args.input, args.output, args.names, args.old_name_field)
    except ValueError as e:
        parser.error(str(e))


if __name__ == '__main__':
    main()
