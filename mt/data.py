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
-- 22.7.2020
   Created Dataset class, replacing code in Aligner.
"""

from .word import *
import os, re, numpy

#from .align import *

DATA_DIR = os.path.join(os.path.dirname(__file__), 'Data')

class Data:
    """
    Data read in from a single file of type .tr
    """

    # regexs to recognize lines in a data file
    languages = re.compile('languages\s*=\s*(.*)')
    comments = re.compile('comments\s*=\s*(.*)')

    # separator for translation groups within a data file
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
        if filename.endswith('.tr'):
            filename = filename[:-3]
        self.filename = filename

    def get_lindex(self, language):
        """
        Get the index of a language the list of languages
        for this Data instance.
        """
        languages = self.info.get('languages')
        if language not in languages:
            print("{} not in dataset languages".format(language))
            return -1
        return languages.index(language)

    def get_words1(self, tgroup, lindices):
        """
        Get the words in a single translation group, given
        language indices lindices.
        """
        return [tgroup.get('words')[lindex] for lindex in lindices]

    def get_words(self, languages):
        """
        Get word sets from the data for languages.
        """
        lindices = [self.get_lindex(language) for language in languages]
        return [self.get_words1(tg, lindices) for tg in self.tgroups]

    @staticmethod
    def create_dataset(datafiles, languages, shuffle=True,
                       test=False, validate=False):
        """
        Create a Data instance for each of datafiles.
        Get the words in all Data instances and create a single
        Dataset instance from these. If shuffle is True,
        randomize the order of the words.
        Create the associated index dataset.
        If test and/or validate is True, split the dataset
        into separate training, test, and validation datasets,
        each with its own index Dataset, returning these.
        """
        wordsets = []
        for datafile in datafiles:
            d = Data(datafile)
            wordsets.extend(d.get_words(languages))
        # base dataset
        ds = Dataset(wordsets, languages, shuffle=shuffle,
                     indices=False, test=test, validate=validate)
        return ds

    # @staticmethod
    # def shuffle(dataset, seedgen=None):
    #     """
    #     Shuffle a dataset so that it can be split into subsets.
    #     If seedgen is non-null,
    #     use seedgen to create the same sets each time this is called.
    #     """
    #     if seedgen:
    #         random.seed(seedgen)
    #     random.shuffle(dataset)

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

        except IOError:
            print('No such data file as {}'.format(path))

    @staticmethod
    def geezify(filename):
        """
        Convert a filename with romanized words to one with geez.
        filename is missing the .tr extension.
        """
        inpath = os.path.join(DATA_DIR, filename + ".tr")
        outpath = os.path.join(DATA_DIR, filename + "ግ.tr")
        with open(inpath, encoding='utf8') as infile:
            with open(outpath, 'w', encoding='utf8') as outfile:
                contents = infile.read().split(Data.tgroup_sep)
                preamble = contents[0]
                tgroups = contents[1:]

                # Write preamble to outfile
                print(preamble, file=outfile)

                # Read and write data
                for tgroup in tgroups:
                    tgroup = tgroup.split('\n')
                    # Each tgroup starts with a comment line
                    print("## {}".format(tgroup[0].strip()), file=outfile)
                    for tg in tgroup[1:]:
                        if not tg:
                            continue
                        rom, cmt = tg.split('#')
                        rootcls, geez = cmt.split(';')
                        cmt = rootcls.strip() + ';' + rom.strip()
                        print("{} # {}".format(geez.strip(), cmt), file=outfile)
#                    if segment:
#                        tgd = [Data.segment(word, Data.ES) for word in tgd]

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
        char0 = word[0]
        # Don't try to segment the word if it's not in the units list
        if char0 not in units[0] and char0 not in units[1]:
            return list(word)
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

class Dataset(list):
    """
    A list of lists of words, one for each of a list of
    languages. Each word is a list of characters or indices
    into character lists.
    """

    # arbitrary constant to generate seed so that datasets are constant
    seedgen = 4
    # proportion of data to use for test set
    test_frac = 0.1
    # proprtion of data to use for validation set
    validation_frac = 0.1
    # constants for different datasets
    TRAINING = 0
    TEST = 2
    VALIDATION = 1

    def __init__(self, words, languages, cat=-1, shuffle=True,
                 parent=None, chars=None, indices=False,
                 test=False, validate=False,
                 constraints=None, searchers=None):
        list.__init__(self, words)
        # a list of language abbreviations, such as amG
        self.languages = languages
        # TRAINING, TEST, or VALIDATION
        self.cat = cat
        # shuffle dataset
        if shuffle:
            self.shuffle(True)
        # parent dataset if this is training, test, or
        # validation set
        self.parent = parent
        if indices:
            # a dataset of character indices, so don't create
            # separate indices dataset
            self.chars = []
            self.indices = None
            self.subds = []
        else:
            # the dataset contains indices rather than characters
            # create the associated index dataset
            self.chars = chars or self.get_all_chars()
            i = self.create_indices()
            self.indices = \
               Dataset(i, languages, cat=cat, shuffle=False,
                       indices=True, chars=None, parent=parent)
            if test or validate:
                subds = self.split(test=test, validate=validate)
                self.subds = subds
            else:
                self.subds = []
        # these are read in from a file with read_constraints()
        self.constraints = constraints
        # these are created in Aligner.make_searchers()
        self.searchers = searchers
        # these are the goal states of searchers after running
        self.alignments = alignments

    def __repr__(self):
        tp = "c" if self.chars else "i"
        lgs = "_".join(self.languages)
        if self.cat >= 0:
            return "{};{};{}".format(tp, self.cat, lgs)
        else:
            return "{};{}".format(tp, lgs)

    def shuffle(self, reproduce=True):
        """
        Shuffle items in dataset so they can be split into training, test,
        validation sets.
        """
        if reproduce:
            numpy.random.seed(Dataset.seedgen)
        numpy.random.shuffle(self)

    def get_sub(self, cat):
        """
        Get the sub-dataset of the given category.
        """
        if self.subds:
            return self.subds[cat]
        return None

    def get_sub_indices(self, cat):
        """
        Get the indices dataset of the sub-dataset
        of the given category.
        """
        sub = self.get_sub(cat)
        if sub:
            return sub.indices
        return None

    def get_all_chars(self):
        """
        Get all of the characters found in the data for each
        language.
        """
        charsets = [set() for l in self.languages]
        for wordset in self:
            for i, word in enumerate(wordset):
                charsets[i].update(word)
        charlists = [list(s) for s in charsets]
        for cl in charlists:
            cl.sort()
        return charlists

    def create_indices(self):
        """
        Create a dataset of indices corresponding to the characters
        in this dataset.
        """
        indices = []
        for wordset in self:
            indices1 = []
            for i, word in enumerate(wordset):
                c_indices = [self.chars[i].index(c) for c in word]
                indices1.append(c_indices)
            indices.append(indices1)
        return indices

    def split(self, test=True, validate=True):
        """Split dataset into test, validation, and training sets."""
        n = len(self)
        test_i = 0
        validation_i = 0
        if test:
            t_i = round(n * Dataset.test_frac)
        if validate:
            v_i = t_i + round(n * Dataset.validation_frac)
        data = self[:t_i], self[t_i:v_i], self[v_i:]
        testset = self.create_subds(data[0], Dataset.TEST)
        valset = self.create_subds(data[1], Dataset.VALIDATION)
        trainset = self.create_subds(data[2], Dataset.TRAINING)
        return trainset, valset, testset

    def create_subds(self, data, cat):
        """
        Create a new sub-dataset from this one with a particular
        category (TRAINING, TEST, VALIDATION).
        """
        return Dataset(data, self.languages, cat=cat,
                       shuffle=False, chars=self.chars,
                       indices=False, parent=self)

    def write(self, filename=None):
        """
        Write segmented words from the dataset to a file.
        (This can then be edited to indicate alignments.)
        """
        filename = filename or self.__repr__() + ".tr"
        path = os.path.join(DATA_DIR, filename)
        with open(path, 'w', encoding='utf8') as file:
            for words in self:
                print(Data.tgroup_sep, file=file)
                for word in words:
                    print("{}".format(' '.join(word)), file=file)

    @staticmethod
    def combine_chars(string):
        """
        Given a string of spaces and characters, return a list
        that combines successive spaces and successive non-spaces.
        """
        stringlist = list(string)
        segs = []
        last = ''
        for char in stringlist:
            if char == ' ':
                if last and ' ' in last:
                    last += ' '
                elif last:
                    segs.append(last)
                    last = ' '
                else:
                    last = ' '
            elif last and ' ' in last:
                segs.append(last)
                last = char
            elif last:
                last += char
            else:
                last = char
        segs += last
        return segs

    @staticmethod
    def wstring2positions(string):
        """
        Given a string of spaces and characters representing
        a word, group the string
        into space and character segments, and return a
        position-keyed dict with within-word positions.
        """
        segs = Dataset.combine_chars(string)
        pos = 0
        chars = {}
        index = 0
        for elem in segs:
            if ' ' not in elem:
                chars[pos] = index
                index += 1
            pos += len(elem)
        return chars

    @staticmethod
    def cstring2positions(string):
        """
        Given a string of spaces and characters representing
        constraints, group the string
        into space and character segments, and return a
        a list of position, constraint element pairs.
        """
        segs = Dataset.combine_chars(string)
        pos = 0
        cons = []
        for elem in segs:
            if ' ' not in elem:
                cons.append((pos, elem))
            pos += len(elem)
        return cons

    @staticmethod
    def record_constraints(spos, tpos, cpos, verbosity=1):
        abs_c = []
        int_c = []
        i_start = None
        if verbosity:
            print("cpos {}".format(cpos))
        positions = set(spos.keys())
        positions.update(set(tpos.keys()))
        positions = list(positions)
        positions.sort()
        current_s = 0
        current_t = 0
        pos_dict = {}
        for p in positions:
            if p in spos:
                current_s = spos.get(p)
            if p in tpos:
                current_t = tpos.get(p)
            pos_dict[p] = (current_s, current_t)
        if verbosity:
            print("pos dict {}".format(pos_dict))
        for pos, label in cpos:
            indices = pos_dict[pos]
            if label == '*':
                # strict positional constraint
                abs_c.append(indices)
            elif label == '[':
                # beginning of an interval
                if pos not in spos:
                    return
                if pos not in tpos:
                    return
                i_start = indices
            elif label == ']':
                # end of an interval; may not correspond to
                # char in source or target
                int_c.append((i_start, indices))
        return abs_c, int_c

    def read_constraints(self, filename=None, verbosity=0):
        """
        Read constraints from a file of aligned words,
        one for each set of words in the dataset.
        """
        filename = filename or self.__repr__() + "A.tr"
        path = os.path.join(DATA_DIR, filename)
        constraints = []
        with open(path, encoding='utf8') as file:
            # skip the empty tgroup before the first ##
            tgroups = file.read().split(Data.tgroup_sep)[1:]
            for tgi, tgroup in enumerate(tgroups):
                sword, tword = self[tgi]
                if verbosity:
                    print("{} ->{}".format(''.join(sword), ''.join(tword)))
                words = tgroup.strip().split("\n")
                # Assume there are three of these:
                # source, target, constraints
                s, t, c = words
                # initialize constraint list
                s_pos = Dataset.wstring2positions(s)
                t_pos = Dataset.wstring2positions(t)
                c_pos = Dataset.cstring2positions(c)
                cons = Dataset.record_constraints(s_pos, t_pos, c_pos,
                                                  verbosity=verbosity)
                if not cons:
                    print("CONSTRAINTS FAILED FOR {}->{}".format(sword, tword))
                    abs_c = int_c = None
                else:
                    abs_c, int_c = cons
                if verbosity:
                    print(" abs constraints: {}".format(abs_c))
                    print(" int constraints: {}".format(int_c))
                constraints.append((abs_c, int_c))
        self.constraints = constraints
#                positions = Dataset.segs2positions(segs, constraint)
#                print(" positions {}".format(positions))
