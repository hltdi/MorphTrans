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
from .utils import *
import os, re, numpy, copy

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
    ES = [["a", "e", "E", "i", "I", "o", "u", "@", "A", "'", "`", ":", "_"],
          {"b": ["b", "bW", "b_", "bW_"],
           "c": ["c", "cW", "c_", "cW_"],
           "C": ["C", "CW", "C_", "CW_"],
           "d": ["d", "dW", "d_", "dW_"],
           "f": ["f", "fW", "f_", "fW_"],
           "g": ["g", "gW", "gY", "g_", "gW_"],
           "h": ["h", "hW", "hY"], "H": ["H", "HW", "HY"],
           "j": ["j", "jW", "j_", "jW_"],
           "k": ["k", "kW", "kY", "k_", "kW_"],
           "K": ["K", "KW", "KY"],
           "l": ["l", "lW", "l_", "lW_"],
           "m": ["m", "mW", "m_", "mW_"],
           "n": ["n", "nW", "n_", "nW_"],
           "p": ["p", "pW", "p_", "pW_"],
           "P": ["P", "PW", "P_", "PW_"],
           "N": ["N", "NW", "N_", "NW_"],
           "q": ["q", "qW", "qY", "q_", "qW_"], "Q": ["Q", "QW"],
           "r": ["r", "rW", "r_", "rW_"],
           "s": ["s", "sW", "s_", "sW_"],
           "S": ["S", "SW", "S_", "SW_"],
           "t": ["t", "tW", "t_", "tW_"],
           "T": ["T", "TW", "T_", "TW_"],
           "w": ["w", "w_"], "y": ["y", "y_"],
           "v": ["v", "vW", "v_"],
           "x": ["x", "xW", "x_", "xW_"],
           "z": ["z", "zW", "z_", "zW_"],
           "Z": ["Z", "ZW", "Z_", "ZW_"],
           "^": ["^s", "^S", "^s_", "^S_", "^h", "^hW",\
                 "^sW", "^SW", "^sW_", "^SW_"]}]

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
        dsid = get_file_id(datafiles[0])
        for datafile in datafiles:
            d = Data(datafile)
            wordsets.extend(d.get_words(languages))
        # base dataset
        ds = Dataset(wordsets, languages, shuffle=shuffle, id=dsid,
                     indices=False, test=test, validate=validate)
        return ds

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
        Convert a file with romanized words to one with geez.
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
    def geminate(filename):
        """
        Convert a file with gemination indicated by ':' to one with
        repeated consonants.
        """
        inpath = os.path.join(DATA_DIR, filename + ".tr")
        outpath = os.path.join(DATA_DIR, filename + "G.tr")
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
                        rom = Data.geminate_word(rom.strip())
                        print("{} # {}".format(rom, cmt.strip()), file=outfile)

    def geminate_dataset(filename):
        """
        Convert a dataset file with gemination indicated by ':' to one with
        repeated consonants.
        """
        inpath = os.path.join(DATA_DIR, filename + ".tr")
        outpath = os.path.join(DATA_DIR, filename + "G.tr")
        with open(inpath, encoding='utf8') as infile:
            with open(outpath, 'w', encoding='utf8') as outfile:
                for line in infile:
                    if '#' in line or '*' in line or '[' in line:
                        print(line, file=outfile, end='')
                    else:
                        # a line with characters
                        word = line.replace('\n', '')
                        word = Data.geminate_word(word)
                        print(word, file=outfile)

    @staticmethod
    def geminate_word(string):
        """
        Replace gemination character with consonant.
        """
        if ':' in string:
            chars = []
            last = ''
            for char in string:
                if char == ':':
                    chars.append(last)
                else:
                    chars.append(char)
                    if char in ["W", "Y"]:
                        last += char
                    elif char != ' ':
                        last = char
            return ''.join(chars)
        return string

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
    # constants for different dataset categories
    TRAINING = 0
    VALIDATION = 1
    TEST = 2
    IDs = ['l', 'v', 't']

    def __init__(self, words, languages, cat=-1, id='',
                 shuffle=True,
                 parent=None, chars=None, indices=False,
                 test=False, validate=False,
                 # these are non-None only for reversed datasets
                 subds=None, nchars=None, index_ds=None,
                 constraints=None, searchers=None, alignments=None):
        list.__init__(self, words)
        # a list of language abbreviations, such as amG
        self.languages = languages
        # TRAINING, TEST, or VALIDATION
        self.cat = cat
        self.id = id
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
            if index_ds:
                self.indices = index_ds
            else:
                i = self.create_indices()
                self.indices = \
                Dataset(i, languages, cat=cat, id=id, shuffle=False,
                        indices=True, chars=None, parent=parent)
            if subds:
                # subdatasets provided only when reversing
                # existing dataset
                self.subds = subds
            elif test or validate:
                subds = self.split(test=test, validate=validate)
                self.subds = subds
            else:
                self.subds = []
        # List of dicts, one for each language.
        # keys: chars, values: count in database.
        # Provided in constructor only when reversing
        # existing database
        self.nchars = nchars or None
        # these are read in from a file with read_constraints()
        # or reversed from original dataset
        self.constraints = constraints
        # these are created in Aligner.make_searchers()
        self.searchers = searchers
        # these are the goal states of searchers after running
        self.alignments = alignments
        # unions of alignments in this dataset and its reverse
        self.union_alignments = None

    ### Names

    def __repr__(self):
        tp = "c" if self.chars else "i"
        lgs = "".join(self.languages)
        catchar = Dataset.catchar(self.cat)
        return "{};{};{};{}".format(self.id, tp, catchar, lgs)

    def merge(self, dataset):
        """
        Add the elements in another dataset to this dataset.
        """
        self.extend(dataset)
        new_chars = dataset.get_all_chars()
        current = self.chars
        chars = []
        for old, new in zip(current, new_chars):
            old = set(old)
            merged = old.union(new)
            merged = list(merged)
            merged.sort()
            chars.append(merged)
        self.chars = chars
        new_indices = dataset.create_indices()
        self.indices.extend(new_indices)

    @staticmethod
    def catchar(catid):
        """
        A character representing the Dataset category.
        """
        if catid < 0:
            return 'd'
        else:
            return Dataset.IDs[catid]

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

    def set_nchars(self, sort=False, add_ds=None):
        """
        Get the numbers of each source and target character
        in the dataset.
        """
        if self.nchars:
            return
        print("Setting nchars dicts for {}".format(self))
        chardicts = [{} for l in self.languages]
        for wordset in self:
            for i, word in enumerate(wordset):
                chardict = chardicts[i]
                for char in word:
                    if char in chardict:
                        chardict[char] += 1
                    else:
                        chardict[char] = 1
        if add_ds:
            # Count characters in additional dataset
            for wordset in add_ds:
                for i, word in enumerate(wordset):
                    chardict = chardicts[i]
                    for char in word:
                        if char in chardict:
                            chardict[char] += 1
                        else:
                            chardict[char] = 1
        if sort:
            for i, cd in enumerate(chardicts):
                cd = list(cd.items())
                cd.sort(key=lambda c: c[1], reverse=True)
                chardicts[i] = cd
        self.nchars = chardicts

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
        testset = self.create_subds(data[0], Dataset.TEST, cons=test)
        valset = self.create_subds(data[1], Dataset.VALIDATION,
                                   cons=validate)
        trainset = self.create_subds(data[2], Dataset.TRAINING)
        return trainset, valset, testset

    def create_subds(self, data, cat, cons=False):
        """
        Create a new sub-dataset from this one with a particular
        category (TRAINING, TEST, VALIDATION).
        """
        dataset = Dataset(data, self.languages, cat=cat, id=self.id,
                          shuffle=False, chars=self.chars,
                          indices=False, parent=self)
        if cons:
            # attempt to read in the constraints for this dataset
            dataset.read_constraints()
        return dataset

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

    ### Reading and creating constraints

    @staticmethod
    def combine_chars(string, verbosity=0):
        """
        Given a string of spaces and characters, return a list
        that combines successive spaces and successive non-spaces.
        """
        stringlist = list(string)
        segs = []
        last = ''
        for char in stringlist:
#            if verbosity:
#                print("char {}, last {}, segs {}".format(char, last, segs))
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
        segs.append(last)
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
    def create_constraints(spos, tpos, cpos, verbosity=1):
        """
        Given source and target positions and constraint
        positions, create the absolute and interval constraint
        lists.
        """
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
        try:
            with open(path, encoding='utf8') as file:
                print("Reading constraints for {} from {}".format(self, path))
                # skip the empty tgroup before the first ##
                tgroups = file.read().split(Data.tgroup_sep)[1:]
                for tgi, tgroup in enumerate(tgroups):
                    sword, tword = self[tgi]
                    if verbosity:
                        print("{} ->{}".format(''.join(sword), ''.join(tword)))
                    words = tgroup.strip().split("\n")
                    # Assume there are three of these:
                    # source, target, constraints
                    if len(words) != 3:
                        print("Problem with {}".format(words))
                    s, t, c = words
                    # initialize constraint list
                    s_pos = Dataset.wstring2positions(s)
                    t_pos = Dataset.wstring2positions(t)
                    c_pos = Dataset.cstring2positions(c)
                    cons = Dataset.create_constraints(s_pos, t_pos, c_pos,
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
        except IOError:
#            print("No constraint file for {}".format(self))
            pass

    ### Reversing a dataset
    def reverse(self, slan_index=1):
        """
        Create a dataset with the first language replaced by
        the language with index slan_index.
        """
        print("Reversing {}".format(self))
        lg = self.languages[:]
        lg[0], lg[slan_index] = lg[slan_index], lg[0]
        id = self.id + 'r'
        items = self.reverse_items(slan_index=slan_index)
        if self.indices:
            # This is a char dataset
            # First create the reversed indices dataset
            index_ds = self.indices.reverse(slan_index=slan_index)
            cons = self.reverse_constraints()
            chars = self.chars[:]
            chars[0], chars[slan_index] = chars[slan_index], chars[0]
            if self.nchars:
                # This may or may not have been set
                nchars = self.nchars[:]
                nchars[0], nchars[slan_index] = nchars[slan_index], nchars[0]
            else:
                nchars = None
            # Create the subdatasets if any
            if self.subds:
                train, val, test = self.subds
                r_train = train.reverse(slan_index)
                r_val = val.reverse(slan_index)
                r_test = test.reverse(slan_index)
                subds = [r_train, r_val, r_test]
            else:
                subds=None
            dataset = \
            Dataset(items, lg, cat=self.cat, id=id,
                    shuffle=False, parent=self.parent, indices=False,
                    chars=chars, nchars=nchars, subds=subds,
                    index_ds=index_ds, constraints=cons)
        else:
            # This is an indices dataset
            dataset = \
            Dataset(items, lg, cat=self.cat, id=id,
                    shuffle=False, parent=self.parent, indices=True)
        # create subds if any
        return dataset

    def reverse_items(self, slan_index=1):
        """
        Create a top-level indices dataset with languages reversed.
        """
        items = []
        for l in self:
            # a list of lists of indices or chars, one for each language
            newl = copy.deepcopy(l)
            newl[0], newl[slan_index] = newl[slan_index], newl[0]
            items.append(newl)
        return items

    def reverse_constraints(self):
        """
        Return a list of constraints, with source and target
        indices swapped.
        """
        if self.constraints:
            print("Reversing constraints for {}".format(self))
            rev_c = []
            for cons in self.constraints:
                abs_c, int_c = cons
                rev_a = []
                rev_i = []
                for sp, tp in abs_c:
                    rev_a.append((tp, sp))
                for beginning, end in int_c:
                    sp0, tp0 = beginning
                    sp1, tp1 = end
                    rev_beg = tp0, sp0
                    rev_end = tp1, sp1
                    rev_i.append((rev_beg, rev_end))
                rev_c.append((rev_a, rev_i))
            return rev_c
        # print("No constraints to reverse!")
        return None

    # character overlap
    def char_overlap(self):
        """
        The proportion of characters in language1 that
        overlap with those in language2.
        """
        overlap = 0
        total = 0
        for wordset in self:
            word0 = wordset[0]
            wset0 = set(word0)
            word1 = wordset[1]
            wset1 = set(word1)
            o1 = wset0.intersection(wset1)
            if o1:
                for char in o1:
                    count = min([word0.count(char), word1.count(char)])
                    overlap += count
            total += len(word0)
        return overlap / total
