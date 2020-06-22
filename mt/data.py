"""
This file is part of the MorphTrans project: https://github.com/hltdi/MorphTrans

    Copyleft 2020. MorphTrans Collaborative.

    MorphTrans is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MorphTrans is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MorphTrans.  If not, see <http://www.gnu.org/licenses/>.

Author: Michael Gasser <gasser@indiana.edu>

-- 30.1.2020
   Created.
"""

import os, re

from .align import *

DATA_DIR = os.path.join(os.path.dirname(__file__), 'Data/CRF')

class Data:

    languages = re.compile('languages\s*=\s*(.*)')
    comments = re.compile('comments\s*=\s*(.*)')

    tgroup_sep = "##"


    # segmentation units for ES languages, used in segment()
    # ES[0] is a list of phones consisting of single characters
    # ES[1] is a dict for characters that may be combined (or not) to make phones
    ES = [["a", "e", "E", "i", "I", "o", "u", "@", "A", "w", "y", "'", "`", ":"],
          {"b": ["b", "bW"], "c": ["c", "cW"], "C": ["C", "CW"],
           "d": ["d", "dW"], "f": ["f", "fW"], "g": ["g", "gW", "gY"],
           "h": ["h", "hW", "hY"], "H": ["H", "HW", "HY"], "j": ["j", "jW"],
           "k": ["k", "kW", "kY"], "K": ["K", "KW", "KY"], "l": ["l", "lW"],
           "m": ["m", "mW"], "n": ["n", "nW"], "p": ["p", "pW"], "P": ["P", "PW"],
           "N": ["N", "NW"], "q": ["q", "qW", "qY"], "Q": ["Q", "QW"],
           "r": ["r", "rW"], "s": ["s", "sW"], "S": ["S", "SW"],
           "t": ["t", "tW"], "T": ["T", "TW"], "v": ["v", "vW"],
           "x": ["x", "xW"], "z": ["z", "zW"], "Z": ["Z", "ZW"],
           "^": ["^s", "^S", "^h", "^hW", "^sW", "^SW"]}]

    def __init__(self, filename, problem=None):
        self.tgroups, self.info = Data.read(filename, problem, segment=True)

    def get_lindex(self, language):
        languages = self.info.get('languages')
        if language not in languages:
            print("{} not in dataset languages")
            return -1
        return languages.index(language)

    def get_words1(self, tgroup, lindices):
        return [tgroup.get('words')[lindex] for lindex in lindices]

    def get_words(self, languages):
        lindices = [self.get_lindex(language) for language in languages]
        return [self.get_words1(tg, lindices) for tg in self.tgroups]

    # New data format: 2020.6.16
    @staticmethod
    def read(filename, problem, segment=True):
        data = []
        info = {}
        path = os.path.join(DATA_DIR, filename)
        try:
            with open(path, encoding='utf8') as datafile:
                contents = datafile.read().split(Data.tgroup_sep)
                preamble = contents[0]
                tgroups = contents[1:]

                # Read preamble
                for line in preamble.split('\n'):
                    line = line.strip()
                    if not line: continue
                    m = Data.languages.match(line)
                    if m:
                        languages = m.group(1).split(',')
                        info['languages'] = languages
                        continue
                    m = Data.comments.match(line)
                    if m:
                        comments = m.group(1).strip()
                        info['comments'] = comments
                        continue
                    print("መስመሩን ማንበብ አልተቻለም።")

                # Read data
                for tgroup in tgroups:
                    tgroup = tgroup.split('\n')
                    tg_data = {}
                    # Each tgroup starts with a comment line
                    tg_data['comments'] = tgroup[0].strip()
                    tgd = []
                    tgc = []
                    for tg in tgroup[1:]:
                        if not tg:
                            continue
                        dt, cmt = tg.split('#')
                        tgd.append(dt.strip())
                        tgc.append(cmt.strip())
                    if segment:
                        tgd = [Data.segment(word, Data.ES) for word in tgd]
                    tg_data['words'] = tgd
                    tg_data['l_comments'] = tgc
                    data.append(tg_data)

            return data, info

#                for line in datafile:
#                    line = line.split('#')[0].strip() # strip comments
#
#                    if not line: continue
#
#                    m = Data.languages.match(line)
#                    if m:
#                        languages = m.group(1).split(',')
#                        info['languages'] = languages
#                        continue
#                    word_list = [word.strip() for word in line.split(',')]
#                    group = Data.create_tgroup(word_list, info.get('languages'), problem)
#                    data.append(group)
#            return data, info
        except IOError:
            print('No such data file as {}'.format(path))

    @staticmethod
    def create_char_list():
        """Make a list of all ES 'characters'."""
        chars = Data.ES[0][:]
        for multichars in Data.ES[1].values():
            chars.extend(multichars)
        return chars

    @staticmethod
    def create_tgroup(lst, languages, problem):
        """lst is a list of strings, one for each word.
        Create an instance of TGroup."""
        group = TGroup([Word(Data.segment(word)) for word in lst], languages, problem)
        return group

    @staticmethod
    def segment(word, units=None, correct=True):
        """Given a string representing a word, segment it into
        a list of phones. Note: list(word) won't work because some phones are
        represented by two characters."""
        units = units or Data.ES
        res = []
        pos = 0
        while pos < len(word):
            ch = word[pos]
            if ch in units[0]:
                res.append(ch)
                pos += 1
            else:
                sublists = units[1]
                if ch in sublists:
                    sublist = sublists[ch]
                    if pos < len(word) - 1 and ch + word[pos + 1] in sublist:
                        if pos < len(word) - 2 and ch + word[pos + 1:pos + 3] in sublist:
                            res.append(ch + word[pos + 1:pos + 3])
                            pos += 3
                        else:
                            res.append(ch + word[pos + 1])
                            pos += 2
                    else:
                        res.append(ch)
                        pos += 1
                elif correct:
                    print(ch, 'in', word, 'is not an acceptable character')
                    return
                else:
                    res.append(ch)
                    pos += 1
        return res
