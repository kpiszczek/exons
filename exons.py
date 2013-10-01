#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division, print_function


import argparse
import re
import csv
import urllib
import mmap
import textwrap

from StringIO import StringIO
from collections import namedtuple

import Bio
# from Bio import Entrez
from Bio.Blast.Applications import NcbiblastpCommandline
# from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
# from Bio import SeqIO


# obiekt porównania dwóch sekwencji
Alignment = namedtuple('Alignment', ['query', 'match', 'subject'])

# dane dlae jedngo aminokwasu z pliku .pdb
PDBRow = namedtuple('PDBRow', ['amino', 'pos', 'struct'])


class __TwoWayDict(dict):

    """
    Słownik dwukierunkowy
    """

    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def update(self, *args, **kwargs):
        if args:
            if len(args) > 1:
                raise TypeError(
                    "update expected at most 1 arguments, got %d" % len(args))
            other = dict(args[0])
            for key in other:
                self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def __len__(self):
        return (dict.__len__(self) + 0.5) // 2

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        dict.__setitem__(self, value, key)

    def __getitem__(self, key):
        if key in self:
            return dict.__getitem__(self, key)
        else:
            return 'X' if len(str(key)) > 2 else 'XXX'


def compareSequences(seq1, seq2):
    """
    compareSequences(seq1, seq2) -> Alignment

    Tworzy obiekt prównania dwóch sekwencji zwartych
    w plikach .fasta.
    """
    output = ''

    try:
        output = NcbiblastpCommandline(
            query=seq1, subject=seq2, outfmt=5, use_sw_tback=True)()[0]
    except Bio.Application.ApplicationError as err:
        print('Brak programu Blast w ścieżce systemowej')
        print(err)
        # print('Próba połączenia się z wersją online...')
        # try:
        #     seq1 = open(seq1).read()
        #     output = NCBIWWW.qblast(
        #         "blastp", 'nt', sequence=seq1, query_file=open(seq2)).read()
        # except IOError:
        #     raise ValueError('Nie udało się pobrać danych')
        exit()

    if not output:
        return

    blast_result_record = NCBIXML.read(StringIO(output))

    alignment = blast_result_record.alignments[0]
    hsp = alignment.hsps[0]
    return Alignment(query=hsp.query, match=hsp.match, subject=hsp.sbjct)


def fetchPDB(pdb_id):
    """
    fetchPDB(pdb_id) -> pdb_file

    Pobiera plik .pdb o zadanym identyfikatorze
    z bazy www.rcsb.org.
    """
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % pdb_id.split('.')[0]
    return urllib.urlopen(url).read()


def getChains(filename):
    try:
        with open(filename) as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

            chains = re.search(
                'CHAIN:\s*(.+);', mm).group(1)

            return chains.replace(' ', '').split(',')
    except IOError:
        print('Nie znaleziono pliku o nazwie "%s"' % filename)
        exit()
    except Exception as err:
        print('Błąd podczas przetwarzania pliku "%s"' % filename)
        print(err)
        exit()


def processPDBFile(filename, chain='A'):
    """
    processPDBFile(filename, chain) -> (sequence, two_dim_struck)

    Tworzy krotkę zwierającą jako pierwszy element
    sekwencje aminokwasów w białku oraz jako drugi element
    słownik zawierający informacje o strukturze 2D, do której
    należy aminokwas na danej pozycji.
    """
    sequence = []
    two_dim_struck = {}
    ids = []
    try:
        with open(filename) as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

            curr_num = ''

            for line in iter(mm.readline, ''):
                if line.startswith('HELIX'):
                    _, _, _, _, c_chain, start, _, _, end, _, _ = line.split()
                    if c_chain == chain:
                        two_dim_struck.update({
                            i: 'H' for i in range(int(start), int(end) + 1)})

                if line.startswith('SHEET'):
                    splitted = line.split()
                    c_chain, start, end = splitted[5], splitted[6], splitted[9]
                    if c_chain == chain:
                        two_dim_struck.update({
                            i: 'B' for i in range(int(start), int(end) + 1)})

                if line.startswith('ATOM'):
                    splitted = line.split()
                    amino, c_chain, num = splitted[3], splitted[4], splitted[5]
                    if c_chain == chain:
                        if num != curr_num:
                            curr_num = num
                            ids.append(int(num))
                            sequence.append(gene_code[amino])

    except IOError:
        print('Nie znaleziono pliku o nazwie "%s"' % filename)
        exit()
    except Exception as err:
        print('Błąd podczas przetwarzania pliku "%s"' % filename)
        print(err)
        exit()

    temp = zip(sequence, ids, (two_dim_struck.get(i, 'L') for i in ids))
    return (''.join(sequence),
           (PDBRow(amino, pos, struct) for amino, pos, struct in temp))


def processGBFile(filename):
    """
    processGBFile(filename) -> (sequence, exons)

    Tworzy krotkę zwierającą jako pierwszy element
    sekwencje aminokwasów kodowaną przez mRNA oraz jako
    drugi element słownik informacje o tym, który ekson
    koduje dany aminokwas.
    """
    sequence = ''
    exons = {}
    try:
        with open(filename) as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

            start, end = re.search(
                'CDS\s*(\d+\.\.\d+)', mm
            ).group(1).split('..')
            start, end = int(start), int(end)

            sequence = re.search(
                '/translation="([^"]+)"', mm).group(1)
            sequence = ''.join(sequence.split())

            exs = re.findall(
                'exon\s*(\d+\.\..*\d+)', mm)
            exs = [(int(x), int(y.replace('<', '')))
                   for x, y in (item.split('..') for item in exs)]

            pos = 0
            for i, (s, e) in enumerate(exs, start=1):
                s = s if s > start else start
                e = e if e < end else end
                for k in range(s, e, 3):
                    pos += 1
                    exons[pos] = i

    except IOError:
        print('Nie znaleziono pliku o nazwie "%s"' % filename)
        exit()
    except Exception as err:
        print('Błąd podczas przetwarzania pliku "%s"' % filename)
        print(err)
        exit()

    return (sequence, exons)


def __makeGeneCodebook():
    """
    __makeGeneCodebook() -> gene_codebook

    Zwaraca słownik dwukierunkowy, tłumaczący kody
    trójliterowe na jednoliterowe oraz odwrotnie.

    >>> gene_codebook = __makeGeneCodebook()
    >>> gene_codebook['GLU'] == 'E'
    >>> gene_codebook['E'] == 'GLU'
    """
    three_letters = (
        'PHE', 'LEU', 'ILE', 'MET', 'VAL', 'SER', 'PRO', 'THR', 'ALA',
        'TYR', 'HIS', 'GLN', 'ASN', 'LYS', 'ASP', 'GLU', 'CYS', 'TRP',
        'ARG', 'GLY')

    one_letter = ('F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
                  'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'G')

    two_way = __TwoWayDict()
    for k, v in zip(three_letters, one_letter):
        two_way[k] = v

    return two_way


gene_code = __makeGeneCodebook()


def __parse_args():
    # parser argmuentów linii komend
    description = '''
    Wizualizajca eksonów w strukturach białek.

    Program do poprawnego działania wymaga dostępu
    do narzędzia BLAST+ z poziomu ścieżki systemowej.
    '''

    epilog = '''
    Program tworzy osobny plik wyjściowy w formacie 'csv',
    dla każdego łańcucha białka zapisanego w pliku .pdb.

    Pliki wyjściowe otrzymują nazwy zgodne ze schematem:
    'output_file_prefix_kod_łańcucha.csv',
    gdzie 'output_file_prefix' jest trzecim parametrem
    wejściowym programu.

    Każdy z plików wyjściowych zaweria pięć kolumn z danymi.
    (id, pdb, match, 2d, exon)
    id - numer kolejnego aminokwasu w dopasowanej sekwencji białka
    pdb - numer aminokwasu w sekwencji z pliku .pdb
    match - kod jednoliterowy aminokwasu
    2d - struktura dwuwymiarowa, w której znajduje się dany aminokwas
    exon - numer eksonu kodującego dany aminokwas
    '''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(description),
        epilog=textwrap.dedent(epilog))

    parser.add_argument('gb_file', help='pełna nazwa pliku .gb')
    parser.add_argument('pdb_file', help='pełna nazwa pliku .pdb')
    parser.add_argument(
        'output_file_prefix', help='prefix nazw plików wyjściowych')

    return parser.parse_args()


def __run():
    """
    Główna funkcja modułu
    """
    args = __parse_args()

    seq2, exons = processGBFile(args.gb_file)
    with open('.gb.fasta', 'w') as f:
        f.write('>%s\n' % (args.gb_file.split('.')[0]))
        f.write('\n'.join(textwrap.wrap(seq2, width=60)))

    for chain in getChains(args.pdb_file):
        seq1, data = processPDBFile(args.pdb_file, chain)
        with open('.pdb.fasta', 'w') as f:
            f.write('>%s\n' % (args.pdb_file.split('.')[0]))
            f.write('\n'.join(textwrap.wrap(seq1, width=60)))

        alignment = compareSequences('.pdb.fasta', '.gb.fasta')

        with open(args.output_file_prefix +
                  '_' + chain + '.csv', 'wb') as csvfile:
            writer = csv.writer(
                csvfile,
                delimiter=',',
                quotechar='|',
                quoting=csv.QUOTE_MINIMAL)

            writer.writerow(['Id', 'PDB', 'matched', '2D', 'exon'])

            match_list = list(alignment.match)
            pdb_data_list = list(data)
            gb_list = list(seq2)
            gb_init = max(
                int(pdb_data_list[0].pos) - 1,
                gb_list.index(pdb_data_list[0].amino))
            i, match_i, pdb_i, gb_i = 0, 0, 0, gb_init
            try:
                while True:
                    amino = match_list[match_i]
                    if amino == ' ':
                        match_i += 1
                        if match_list[match_i] == pdb_data_list[pdb_i
                                                                + 1].amino:
                            pdb_i += 1
                        if match_list[match_i] == gb_list[gb_i + 1]:
                            gb_i += 1
                        else:
                            writer.writerow([i + 1, '-', '-', '-', '-'])
                            i += 1
                            continue
                    else:
                        dat = pdb_data_list[pdb_i]

                    row = ['_'] * 5
                    row[0] = i + 1

                    if amino == dat.amino:
                        row[1], row[2], row[3] = dat.pos, dat.amino, dat.struct
                    if amino == gb_list[gb_i]:
                        row[4] = exons[gb_i + 1]
                    writer.writerow(row)
                    i += 1
                    match_i += 1
                    pdb_i += 1
                    gb_i += 1
            except IndexError:
                pass


if __name__ == "__main__":
    __run()
