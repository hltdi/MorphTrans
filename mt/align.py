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

-- 19.6.2020
   Created.
   Aligner uses the EM algorithm and word translation pairs in two languages
   (for convenience, referred to as 'source' and 'target') to learn alignments
   of characters in the words and character-to-character probabilities.
"""

import os
from .data import *
from .utils import *
from .search import *

class Aligner:

    # probability of the same character in both languages
    eqprob = 0.5
    # other probabilities
    otherprob = 0.01
    # null character
    NULL = '0'
    # number of length differences for parameters
    nlendiffs = 16
    # smallest length difference for parameters
    lendiff0 = -4
    # arbitrary constant to generate seed so that datasets are constant
    seedgen = 4
    # proportion of data to use for test set
    test_frac = 0.1
    # proprtion of data to use for validation set
    validation_frac = 0.1
    # constants for different datasets
    # TRAINING = 0
    # TEST = 1
    # VALIDATION = 2
    # IDs for different datasets
    # DS_IDS = ['t', '?', 'v']

    # Directory for storing parameter dumps
    dump_dir = os.path.join(os.path.dirname(__file__), 'Dumps')

    def __init__(self, source, target, datafiles=None,
                 test=True, validate=True,
                 restore=False):
        self.source = source
        self.target = target
        # Source and target language characters
        self.schars = self.get_schars()
        self.tchars = self.get_tchars()
        # History of sum of squares of changes in probs during M-step
        self.changes = []
        ## Variables depending on data
        # Subsets of data for training, test, validation
        # List of word pairs, consisting of character strings
        self.data = []
        # List of word pairs, consisting of lists of character indices
        self.training = []
        self.test = []
        self.validation = []
        ## Other arrays needed for training
        # Arrays that depend on data
        self.counts = []
        self.tchar_counts = []
        self.tdelt_counts = []
        self.ljustcounts = []
        self.rjustcounts = []
        self.diffcounts = []
        self.tindex_combs = []
        if datafiles:
            self.data = Data.create_dataset(datafiles,
                                            [self.source, self.target],
                                            test=test, validate=validate)
            # get source and target characters in data
            self.schars, self.tchars = self.data.chars
            # Number of source and target characters
            self.ns = len(self.schars)
            self.nt = len(self.tchars)
            ## Trained parameters (saved with save(), loaded with load())
            # Character association probabilities
            self.make_probs()
            # self.probs = self.make_probs()
            # self.init_probs()
            # Probabilities of deleting target characters
            self.tdeltas = []
            # Justification/alignment parameters
            self.lambdas = self.make_ljust_probs()
            # array to store normalization denominator in M-step
            self.tnorms = numpy.zeros(self.nt+1)
            # Shuffle data so that it can be split into training, testing
            # and validation subsets
            # self.shuffle(self.data)
            if test or validate:
                self.training, self.validation, self.test = self.data.subds
            else:
                self.training = self.data
            self.trainingN = len(self.training)
            self.counts = self.make_counts()
            # Target character deletion parameters and counts
            self.tdelt_counts = self.make_tdel_counts()
            # Number of target characters
            self.tchar_counts = self.make_tchar_counts()
            # Probabilities of deleting target characters
            self.tdeltas = self.make_tdeltas()
            # Store combinations of target indices to avoid
            # repeating itertools.combinations()
            self.tindex_combs = self.make_target_index_combs()
            # Alignment/justification counts
            self.ljustcounts = self.make_just_counts()
            self.rjustcounts = self.make_just_counts()
            self.diffcounts = self.make_diff_counts()
            # Identifier for Aligner based on filenames
            self.id = get_file_id(datafiles[0])
        else:
            # no data provided; set parameters on the basis
            # of list of all source and target characters
            self.ns = len(self.schars)
            self.nt = len(self.tchars)
            self.make_probs()
            self.tdeltas = []
            self.lambdas = self.make_ljust_probs()
            self.tnorms = numpy.zeros(self.nt+1)
            self.id = ''

    def __repr__(self):
        return "↙↘:{}:{},{}".format(self.id, self.source, self.target)

    ### experiments

    def prepare(self, validate=False):
        """
        Prepare for an experiment.
        """
        # create constraints
        #   assumes there is a test or validation dataset
        #   and a constraint file
        if validate:
            dataset = self.validate
            dsid = Dataset.VALIDATION
        else:
            dataset = self.test
            dsid = Dataset.TEST
        print("Preparing run")
        if not dataset.constraints:
            dataset.read_constraints()
        # create searchers for each word pair
        if not dataset.searchers:
            self.make_searchers(dsid)

    def run(self, validate=False, iter_cutoff=8, restore=None):
        """
        Run an alignment experiment on the given dataset.
        """
        dscat = Dataset.VALIDATION if validate else Dataset.TEST
        restored = False
        if restore:
            restored = self.restore(restore)
        self.prepare(validate=validate)
        if not restored:
            # train if restore fails too
            self.EM(iter_cutoff=iter_cutoff)
        self.align(dscat)

    def dump(self):
        """
        Dump the trained parameters to files.
        """
        time = time2str(get_time(), short=True, file=True)
        filename = "{}_{}.npy".format(self.id, time)
        print("Dumping parameters to {}".format(filename))
        path = os.path.join(Aligner.dump_dir, filename)
        with open(path, 'wb') as file:
            numpy.save(file, self.probs)
            numpy.save(file, self.tdeltas)
            numpy.save(file, self.lambdas)

    def restore(self, filename):
        """
        Restore trained parameters from a file.
        """
        print("Attempting to restore parameters from {}".format(filename))
        path = os.path.join(Aligner.dump_dir, filename)
        try:
            with open(path, 'rb') as file:
                self.probs = numpy.load(file, allow_pickle=True)
                self.tdeltas = numpy.load(file, allow_pickle=True)
                self.lambdas = numpy.load(file, allow_pickle=True)
            return True
        except IOError:
            print("Unable to read from {}".format(path))
            return False

    ### Datasets and words

    def get_data_indices(self, training=False, test=False, validate=False):
        """
        Word pairs as indices into various parameter arrays.
        """
        if training:
            return self.training.indices
        elif test:
            return self.test.indices
        elif validate:
            return self.validation.indices
        else:
            return self.data.indices

    def get_data_by_cat(self, cat):
        """
        Get the dataset given an integer category.
        """
        if cat == Dataset.TRAINING:
            return self.training
        elif cat == Dataset.TEST:
            return self.test
        elif cat == Dataset.VALIDATION:
            return self.validation
        else:
            return self.data

    def get_words_by_di(self, data_i, dscat=-1, data=None):
        """
        Given an index into a dataset and possibly
        a dataset category, return the word indices,
        word characters, and constraints for that
        index.
        """
        data = data or self.get_data_by_cat(dscat)
        wp_indices = data.indices[data_i]
        wp_chars = data[data_i]
        cons = data.constraints[data_i] if data.constraints else None
        return wp_indices, wp_chars, cons

    ### Accessing characters and character indices

    def get_schars(self, every=False):
        """
        All characters for source language.
        """
        phones = Phone.PHONES.get(self.source)
        if every:
            return phones[0] + phones[1]
        return phones[0]

    def get_tchars(self, every=False):
        """
        All characters for target language.
        """
        phones = Phone.PHONES.get(self.target)
        if every:
            return phones[0] + phones[1]
        return phones[0]

    def get_tchar_index(self, char, check=False):
        """
        Index associated with target character.
        """
        if check:
            if char not in self.tchars:
                return -1
        return self.tchars.index(char)

    def get_schar_index(self, char, check=False):
        """
        Index associated with source character.
        """
        if check:
            if char not in self.schars:
                return -1
        return self.schars.index(char)

    def get_schar(self, index):
        """
        Character associated with source character index.
        """
        if index >= self.ns:
            return Aligner.NULL
        return self.schars[index]

    def get_tchar(self, index):
        """
        Character associated with target character index.
        """
        if index >= self.nt:
            return Aligner.NULL
        return self.tchars[index]

    ### Accessing parameter values

    def get_prob(self, sindex, tindex):
        """
        Connection probability for source character
        with index sindex and target character with
        index tindex.
        """
        return self.probs[sindex][tindex]

    def get_schar_probs(self, char):
        """
        All connection probabilities associated with
        source character char.
        """
        charindex = self.get_schar_index(char)
        return self.probs[charindex]

    def get_tchar_probs(self, char):
        """
        All connection probability associated with
        target character char.
        """
        charindex = self.get_tchar_index(char)
        return self.probs[:,charindex]

    def get_char_prob(self, schar, tchar):
        """
        Connection probability for source character
        schar and target character tchar.
        """
        sindex = self.get_schar_index(schar)
        tindex = self.get_tchar_index(tchar)
        return self.get_prob(sindex, tindex)

    def get_target_probs(self, tword, sindex):
        """
        Given a sindex, find the connection probabilities
        of all characters in tword, with NULL in last position.
        """
        tprobs = [self.probs[sindex,tindex] for tindex in tword]
        tprobs += [self.probs[sindex,self.nt]]
        return numpy.array(tprobs)

    def get_tdeltas(self, tword):
        """
        Return an array with all of the target delta
        probabilities for characters in tword.
        """
        tdeltas = [self.tdeltas[tindex] for tindex in tword]
        return numpy.array(tdeltas)

    def get_ljust_prob(self, sword, tword):
        """
        Left justification probability for source word
        sword and target word tword (depending on length
        difference).
        """
        m = len(sword)
        l = len(tword)
        index = self.get_lendiff_index(m-l)
        return self.lambdas[index]

    def get_NULL_prob(self, si):
        """
        Probability that source character with index si
        is deleted (not connected to any target character).
        """
        return self.probs[si,self.nt]

    def get_lendiff(self, index):
        """
        The len(s) - len(t) difference associated with the
        count or parameter in the index position.
        """
        return Aligner.lendiff0 + index

    def get_lendiff_index(self, diff):
        """
        The index associated with the len(s) - len(t) difference
        in the count or parameter table.
        """
        return diff - Aligner.lendiff0

    def targetN(self):
        """
        Number of target probs and counts, with the NULL position
        at the end.
        """
        return self.nt+1

    ### Initializing parameters and other arrays.

    def make_probs(self):
        """
        Make the table of character connection probabilities:
        p(s|t) for s/t character combinations.
        source is rows, target columns.
        """
        array = numpy.full((self.ns, self.targetN()), Aligner.otherprob)
        # initialize the array
        getrowchar = self.get_schar
        getcolindex = self.get_tchar_index
        for i, row in enumerate(array):
            char = getrowchar(i)
            colindex = getcolindex(char, True)
            if colindex >= 0:
                array[i,colindex] = Aligner.eqprob
            # Deletion probability
            array[i,self.nt] = 0.2
        self.probs = array
        # Normalize
        self.norm_cols()
#        return array

    # def init_probs(self):
    #     """
    #     Initialize the connection probabilities in the table.
    #     """
    #     array = self.probs

    def norm_cols(self):
        """
        Normalize each column of the probs array except the last
        one (deletion).
        """
        array = self.probs
        for c in range(array.shape[1]-1):
            self.norm_col(c, array)

    def norm_col(self, ci, array):
        """
        Normalize probabilities in column ci so that they sum to 1.0.
        """
        # sum of pross in column ci
        colsum = array[:,ci].sum()
        # divide each prob in c1 by colsum
        if colsum == 0.0:
            print("Can't normalize 0 column")
        else:
            array[:,ci] /= colsum

    def make_counts(self):
        """
        Make the 'st' table of normalized counts: c(s|t) for s/t character
        connections. source is rows, target columns.
        """
        array = numpy.zeros((self.ns, self.targetN(), self.trainingN))
        return array

    def make_just_counts(self):
        """
        Make the left- or right-justified table of normalized counts,
        one for each combination of len(s) - len(t) distance."""
        array = numpy.zeros((Aligner.nlendiffs))
        return array

    def make_diff_counts(self):
        """
        Count the number of pairs for each len(s) - len(t) difference.
        """
        array = numpy.zeros((Aligner.nlendiffs))
        for wpair in self.get_data_indices(training=True):
            diff = len(wpair[0]) - len(wpair[1])
            index = self.get_lendiff_index(diff)
            array[index] += 1
        return array

    def make_tdeltas(self):
        """
        Make the parameters governing whether characters in a
        target word will be deleted (not connected to any
        source character).
        """
        array = numpy.zeros((self.nt))
        for tindex in range(self.nt):
            if self.tchar_counts[tindex] > 0:
                # we only need values for target characters
                # that actually occur in the dataset
                array[tindex] = 0.1
        return array

    def make_target_index_combs(self):
        """
        Create a set of index combinations for each target word
        length.
        """
        tlengths = {len(wp[1]) for wp in self.get_data_indices(training=True)}
        combs = [(n, all_index_combs(n)) for n in tlengths]
        return dict(combs)

    def make_tchar_counts(self):
        """
        Count the number of times each target character
        occurs in the dataset.
        """
        occs = []
        for index in range(self.nt):
            occs.append(sum([wp[1].count(index)\
                        for wp in self.get_data_indices(training=True)]))
        return numpy.array(occs)

    def make_tdel_counts(self):
        """
        Make the array to store weighted counts of deletions
        for each target character.
        """
        array = numpy.zeros((self.nt))
        return array

    def make_ljust_probs(self):
        """
        Make the probabilities of left alignment, given different
        differences between source and target word lengths.
        """
        array = numpy.full((Aligner.nlendiffs), 0.5)
        return array

    def init_lambdas(self):
        """
        Initialize the lambdas (left justification parameters).
        """
        self.lambdas[:] = 0.5

    ### EM

    def EM(self, iter_cutoff=10, error_cutoff=-20.0, verbosity=0):
        """
        Run the EM algorithm for up to iter_cutoff iterations
        or up to a RMS error of error_cutoff.
        """
        iteration = 0
        changes = 1.0
        while iteration < iter_cutoff and changes > error_cutoff:
            print("ITERATION {}".format(iteration+1))
            print("E step")
            self.e_init()
            self.update_counts(verbosity=verbosity)
            print("M step")
            self.update_tnorms(verbosity=verbosity)
            self.update_probs(verbosity=verbosity)
            self.update_lambdas(verbosity=verbosity)
            self.update_tdeltas(verbosity=verbosity)
            iteration += 1
            changes = self.changes[-1]

    def reinit(self):
        """
        Reinitialize parameters for a new EM run.
        """
        self.init_probs()
        self.init_lambdas()

    ## E-step

    def e_init(self):
        """
        Do whatever initialization is required before the E-step
        on each EM iteration.
        """
        # Reset ljust and rjust counts to 0
        self.ljustcounts[:] = 0.0
        self.rjustcounts[:] = 0.0
        # Reset target delta counts to 0
        self.tdelt_counts[:] = 0.0

    def update_counts(self, verbosity=0):
        """
        On the E-step, update all of the expected counts for each
        word pair, based on current parameters, estimated
        cooccurrences of chars and estimated deleted characters
        in target and source words.
        """
        for wpindex, wpair in enumerate(self.get_data_indices(training=True)):
            sword, tword = wpair
            slength = len(sword)
            tlength = len(tword)
            lendiff = slength - tlength
            diffindex = self.get_lendiff_index(lendiff)
            if verbosity:
                print("Word {}, m {}, l {}".format(wpindex, slength, tlength))
            # List of target delta parameters for this word
            tdelts = self.get_tdeltas(tword)
            # Estimated number of target char deletions based on
            # tdel probabilities
            tdels = self.calc_exp_tdels(tword, tdelts)
            # Estimated number of target char deletions adjusted
            # for word length difference
            adj_tdels = max([-lendiff, tdels])
            if verbosity:
                print(" exp tdels {}, adjusted {}".format(tdels, adj_tdels))
            sdels = self.calc_exp_sdels(sword, tword, adj_tdels)
            if verbosity:
                print(" exp sdels {}".format(sdels))
            for si in sword:
                denom = self.calc_prob_sum(si, tword)
                # Aligned with characters in target word
                for ti in tword:
                    self.update_count(si, ti, sword, tword, wpindex,
                                      denom, diffindex)
                self.update_sdel_count(si, sword, tword, wpindex,
                                       expdels=sdels, denom=denom)
            self.update_tdelta_counts(tword, tdelts, adj_tdels,
                                      verbosity=verbosity)

    def update_count(self, si, ti, sword, tword, wpindex, denom=None,
                     diffindex=0):
        """
        Calculate the current expected weighted count of the character
        index si associated with the character index ti, given
        ti.
        And increment the left and right justification counts.
        """
        if not denom:
            # Denominator is the s|t probability normalized by the
            # s|Σt_i for all i in target word.
            denom = self.calc_prob_sum(si, tword)
        prob = self.probs[si,ti]
        leftprox, rightprox = self.calc_cooccs(sword, si, tword, ti)
        if leftprox or rightprox:
            cooccs = (leftprox + rightprox) / 2.0
            self.counts[si, ti, wpindex] = (prob * cooccs) / denom
            # also normalize these?
            self.ljustcounts[diffindex] += (prob * leftprox)
            self.rjustcounts[diffindex] += (prob * rightprox)

    def calc_prob_sum(self, si, tword):
        """
        Calculate the sum of all probs with s = si and ti in tword.
        Normalizing factor for updating probabilities.
        """
        probs = [self.probs[si,ti] for ti in tword]
        return sum(probs) + self.get_NULL_prob(si)

    def calc_cooccs(self, sword, schar, tword, tchar):
        """
        Product of number of schars in sword and number of tchars in tword
        weighted by left justification probabilities applied to characters'
        positions.
        Return both left and right estimates.
        """
        total = 0.0
        lbd = self.get_ljust_prob(sword, tword)
        slength = len(sword)
        tlength = len(tword)
        maxsep = max([slength, tlength])
        leftproxsum = 0.0
        rightproxsum = 0.0
        for j, s in enumerate(sword):
            if s != schar:
                continue
            for i, t in enumerate(tword):
                if t != tchar:
                    continue
                # characters at positions j, j match schar, tchars
                sep = j - i
                leftsep = abs(sep)
                rightsep = abs(j-slength - (i - tlength))
                leftsepnorm = leftsep / maxsep
                rightsepnorm = rightsep / maxsep
                leftprox = lbd * (1.0 - leftsepnorm)
                rightprox = (1.0 - lbd) * (1.0 - rightsepnorm)
                total += leftprox + rightprox
                # update justification make_counts
                leftproxsum += leftprox
                rightproxsum += rightprox
        return leftproxsum, rightproxsum

    # Deletion counts on E-step

    def update_sdel_count(self, si, sword, tword, wpindex,
                          expdels=0, denom=None):
        """
        Calculate the current expected weighted count of deletions
        of character si for word pair wpair.
        """
        if not denom:
            denom = self.calc_prob_sum(si, tword)
        prob = self.get_NULL_prob(si)
        m = len(sword)
        l = len(tword)
        # Number of instances of si in sword
        scount = sword.count(si)
        # Proportion of characters in sword that are si
        sprop = scount / m
        # Expected number of deletions of si
        si_dels = sprop * expdels
        #deletions
        if si_dels:
            value = (prob * si_dels)
            self.counts[si, self.nt, wpindex] = value / denom

    def update_tdelta_counts(self, tword, tdelts, exptdels,
                             verbosity=0):
        """
        Update the counts of scaled deletions of target
        characters in the current word.
        """
        if exptdels == 0.0:
            return
        tdeltsum = sum(tdelts)
        for tindex, tdelt in zip(tword, tdelts):
            exptdel = tdelt * exptdels / tdeltsum
            self.tdelt_counts[tindex] += exptdel
            if verbosity > 1:
                print("  Updating tdcount for {}: {}".format(tindex,
                self.tdelt_counts[tindex]))

    def calc_exp_tdels(self, tword, tdelts):
        """
        Calculate the expected number of target char
        deletions, based on their character deletion
        probabilities.
        """
        icombs = self.tindex_combs[len(tword)]
        return weighted_seq_prob_combs(tdelts, icombs)

    def calc_exp_sdels(self, sword, tword, tdels):
        """
        Calculate the expected number of source char
        deletions, based on the expected target char deletions
        and the lengths of the words.
        """
        ldiff = len(sword) - len(tword)
        return max([0, ldiff + tdels])

    ## M-step

    def update_tnorms(self, verbosity=0):
        """
        Using current counts, calculate the normalizing factor
        for each column (λ_e in Brown et al.).
        """
        for tindex in range(self.targetN()):
            total = self.counts[:,tindex,:].sum()
            self.tnorms[tindex] = total

    def update_probs(self, verbosity=0):
        """
        Update all of the probabilities in the M-step.
        """
        # This might not need to be
        changes = 0.0
        for sindex in range(self.ns):
            for tindex in range(self.targetN()):
                oldprob = self.probs[sindex,tindex]
                total = self.counts[sindex,tindex,:].sum()
                if total:
                    self.probs[sindex,tindex] = total / self.tnorms[tindex]
                else:
                    self.probs[sindex,tindex] = 0.0
                change = numpy.power(self.probs[sindex,tindex] - oldprob, 2)
                changes += change
        changes /= self.probs.size
        changes = numpy.log(changes)
        self.changes.append(changes)
        print("Changes: {}".format(changes))

    def update_lambdas(self, verbosity=0):
        """
        Update all left justification probabilities.
        """
        for lendiff in range(Aligner.nlendiffs):
            # Find the left and right justification counts for this lendiff
            leftcount = self.ljustcounts[lendiff]
            rightcount = self.rjustcounts[lendiff]
            if leftcount or rightcount:
                self.lambdas[lendiff] = leftcount / (leftcount + rightcount)

    def update_tdeltas(self, verbosity=0):
        """"
        Update the target character deletion parameters.
        """
        for tindex, td_count in enumerate(self.tdelt_counts):
            if td_count:
                t_sum = self.tchar_counts[tindex]
                norm = td_count / t_sum
                if verbosity > 1:
                    print("  New tdelt for {}: {}".format(tindex, norm))
                self.tdeltas[tindex] = norm

    ### Relative positions within source and target Words

    # @staticmethod
    # def pdf(x, mu=0.0, sigma=0.0):
    #     sigma = sigma or Aligner.sigma
    #     x = float(x - mu) / sigma
    #     return numpy.exp(-x*x/2.0) / numpy.sqrt(2.0*numpy.pi) / sigma

    # @staticmethod
    # def rel_index(index, length):
    #     """
    #     Return the relative position of index within word.
    #     IF index = 0 => 0.0
    #     IF index = len(word)-1 => 1.0
    #     """
    #     # if index < 0:
    #     #     # convert negative index to positive
    #     #     index = length + index
    #     # n = length-1
    #     return index / (length - 1)
    #
    # @staticmethod
    # def rel_sep(length1, index1, length2, index2):
    #     """
    #     The relative left-justified separation of index1 in word of length1
    #     and index2 in word of length2.
    #     0.0 <= rvalue <= 1.0
    #     """
    #     relindex1 = Aligner.rel_index(index1, length1)
    #     relindex2 = Aligner.rel_index(index2, length2)
    #     return abs(relindex1 - relindex2)
    #
    # @staticmethod
    # def right_rel_sep(length1, index1, length2, index2):
    #     """
    #     The relative right-justified separation of index1 in word of length1
    #     and index2 in word of length2.
    #     0.0 <= rvalue <= 1.0
    #     """
    #     relindex1 = Aligner.rel_index(length1 - index1, length1)
    #     relindex2 = Aligner.rel_index(length2 - index2, length2)
    #     return abs(relindex1 - relindex2)
    #
    # @staticmethod
    # def rel_proximity_sum(sindex, sword, tword):
    #     l = len(sword)
    #     m = len(tword)
    #     return sum([Aligner.scaled_rel_proximity(l, sindex, m, tindex) \
    #     for tindex in range(m)])
    #
    # @staticmethod
    # def rel_proximity_prob(sindex, sword, tindex, tword, total=0):
    #     if not total:
    #         total = Aligner.rel_proximity_sum(sindex, sword, tword)
    #     l = len(sword)
    #     m = len(tword)
    #     proximity = Aligner.scaled_rel_proximity(l, sindex, m, tindex)
    #     return proximity/l
    #
    # @staticmethod
    # def rel_proximity(word1, index1, word2, index2, length1=0, length2=0):
    #     """
    #     Relative proximity of position index1 in word1 and index2
    #     in word2, scaled by delta.
    #     """
    #     length1 = length1 or len(word1)
    #     length2 = length2 or len(word2)
    #     distance = \
    #     Aligner.rel_dist(length1, index1, length2, index2)
    #     return 1.0 - Aligner.delta * distance

    # @staticmethod
    # def rel_proximity_sum(sindex, sword, tword):
    #     l = len(sword)
    #     m = len(tword)
    #     return sum([Aligner.rel_proximity(sword, sindex, tword, tindex, l, m) \
    #     for tindex in range(m)])

    ### Alignments for word pairs and for a whole dataset.

    def align(self, dscat=-1):
        """
        Create searchers if they don't already exist.
        Run each searcher and save the resulting goal alignment.
        """
        data = self.get_data_by_cat(dscat)
        print("Aligning {}".format(data))
        if not data.searchers:
            self.make_searchers(dscat=dscat)
        alignments = []
        for searcher in data.searchers:
            result = searcher.run()
            if result:
                goal_alignment, cost = result
                alignments.append(goal_alignment)
            else:
                alignments.append(None)
        data.alignments = alignments

    def make_searchers(self, dscat=-1):
        """
        Create searcher and alignment objects for each word pairs in the
        indicated dataset.
        """
        data = self.get_data_by_cat(dscat)
        print("Creating searchers for {}".format(data))
        dsid = Dataset.IDs[dscat] if dscat >=0 else ''
        searchers = []
        for d_i in range(len(data)):
            wpi, wpc, cons = self.get_words_by_di(d_i, data=data)
#            wp = data.indices[d_i]
#            wp_chars = data[d_i]
#            cons = data.constraints[d_i]
            name = "Aligner_{}{}".format(dsid, d_i)
            searcher = self.make_searcher(name, wpi, wpc, cons)
            searchers.append(searcher)
        data.searchers = searchers

    def make_alignment(self, wpi=None, wpc=None, cons=None,
                       direction=0,
                       data_i=None, explicit=None, dscat=-1):
        """
        Create an Alignment object for the wpair at index data_i.
        """
        if not wpi:
            wpi, wpc, cons = self.get_words_by_di(data_i, dscat)
        return Alignment(wpi, wpc,
                         self, explicit=explicit, constraints=cons,
                         direction=direction)

    def make_searcher(self, name, wpi=None, wpc=None, cons=None,
                      data_i=-1, dscat=-1, direction=0,
                      verbosity=0):
        """
        Create a best-first searcher for alignments, using
        instances of the Alignment class as states.
        """
        # Use the data index and dataset category to get the
        # words
        if not wpi:
            wpi, wpc, cons = self.get_words_by_di(data_i, dscat)
        return \
        BestFirst(name,
                  goal_test=lambda s: s.complete(),
                  make_start=lambda: self.make_alignment(wpi, wpc, cons, direction),
                  extend=lambda s: s.extend(verbosity=verbosity),
                  evaluate=lambda s: s.cost + s.distance)

class Alignment(list):
    """
    Alignment between a source word and its translation in
    one target language, implemented as a subclass of list,
    with indices into positions in the target word as list elements.
    """

    beam = 5
    '''Number of succeeding positions to align with current
    s position when extending.'''

    def __init__(self, wpair, wpair_chars, aligner,
                 cost=0, distance=0, explicit=None,
                 direction=0, constraints=None):
        # explicit is one or more target positions specified
        # at initialization
        if explicit:
            list.__init__(self, explicit)
        else:
            list.__init__(self)
        # source and target words as lists of integers
        self.sword = wpair[0]
        self.tword = wpair[1]
        # lengths of source and target words
        self.slength = len(self.sword)
        self.tlength = len(self.tword)
        # associated Aligner instance; neede because it's got all
        # of the parameters and characters
        self.aligner = aligner
        # source and target words as lists of character strings
        self.schars = wpair_chars[0]
        self.tchars = wpair_chars[1]
        # accumulated evaluation
        self.cost = cost
        # number of source characters left to align
        self.distance = distance or len(self.sword)
        # whether to set alignment direction now
        if direction == 0:
            self.set_direction()
        else:
            self.left = True if direction == 1 else False
        # constraints on valid alignments for this word pair
        self.constraints = constraints

    def __repr__(self):
        #return "{}→{}\n{}".format(
        #''.join(self.schars),
        #''.join(self.tchars),
        return self.pretty()

    def show(self):
        """
        Because __repr__() obscures the raw list representation,
        we need this to show it.
        """
        return tuple(self)

    ### Access

    def swi(self, index):
        """
        The source word index associated with position index.
        If the alignment is right-justified, this represents
        the distance from the right end of the word.
        """
        if not self.left:
            return self.sword[-index-1]
        return self.sword[index]

    def twi(self, index):
        """
        The target word index associated with position index,
        via this alignment.
        """
        if not self.left:
            return self.tword[-index-1]
        return self.tword[index]

    def schar(self, index):
        """
        The source word character associated with position
        index.
        """
        if not self.left:
            return self.schars[-index-1]
        return self.schars[index]

    def tchar(self, index):
        """
        The target word character associated with position index,
        via this alignment.
        """
        if not self.left:
            return self.tchars[-index-1]
        return self.tchars[index]

    def last_nondel(self):
        """
        For left alignment, the rightmost non-deleted aligned position,
        for right alignment, the leftmost non-deleted aligned position
        (though both are searched for from the end of the self list),
        -1 if there are none.
        """
        nondel = [p for p in self if not self.deleted(p)]
        if nondel:
            return nondel[-1]
        else:
            return -1

    ### Alignment direction, based on justification parameters

    def set_direction(self):
        """
        Set the direction of alignment after training.
        """
        lam = self.aligner.get_ljust_prob(self.sword, self.tword)
        if lam >= 0.5:
            self.left = True
        else:
            self.left = False

    def right2left(self):
        """
        Convert an alignment starting on the right to one
        starting on the left.
        """
        leftalg = [-1] * self.slength
        for spos, tpos in enumerate(self):
            if tpos != -1:
                leftalg[-spos-1] = self.tlength-tpos-1
        return self.copy(explicit=leftalg, direction=1)

    ### Extending: creating new search states

    def copy(self, newcost=0, direction=-1,
             explicit=None):
        """
        Return a copy of the Alignment with the same
        word pair and indices, cost, and distance.
        Explicit is set to the current alignment values
        for this alignment.
        """
        if direction < 0:
            direction = 1 if self.left else 2
        return Alignment((self.sword, self.tword),
                         (self.schars, self.tchars),
                         self.aligner,
                         explicit=explicit or self[:],
                         cost=self.cost + newcost,
                         distance=self.distance,
                         direction=direction,
                         constraints=self.constraints)

    def update(self, tpos):
        """
        Update an alignment copy with the latest index.
        """
        self.append(tpos)
        self.update_cost()
        self.distance -= 1

    def update_cost(self):
        """
        Add cost of last index.
        """
        index = len(self) - 1
        prob = self.prob1(index) or 1.0e-20
        self.cost -= numpy.log(prob)

    def deleted(self, position):
        """
        Is the character in this source word position
        deleted (not connected to any target word character?
        """
        return position < 0

    def complete(self):
        """
        Is the alignment complete? Does it have a connection
        or deletion for every character?
        """
        return len(self) == self.slength

    def next(self):
        """
        The character following the last one aligned.
        """
        if self.complete():
            # No more characters
            return
        return self.swi(len(self))

    def extend(self, verbosity=0):
        """
        Create new alignments starting with the current one by
        trying different connections for the next character
        and deletion of the next character.
        This will produce a maximum of Alignment.beam + 1
        new alignments.
        """
        if verbosity:
            print("Extending {}".format(tuple(self)))
        if self.complete():
            if verbosity:
                print(" Complete")
            # impossible to extend a complete state
            return
        alignments = []
        # make new alignments of the current source position
        # with the next beam positions, starting at the
        # last target position connected to a source character
        if len(self) == 0 or all([self.deleted(p) for p in self]):
            last_tpos = -1
        else:
            last_tpos = self.last_nondel()
        if verbosity:
            print(" Last tpos {}".format(last_tpos))
        # start connecting here
        tpos = last_tpos + 1
        # number of target positions tried
        beam = 0
        # cost to add to new alignment because of deleted target character
        newcost = 0
        while tpos < self.tlength and beam < Alignment.beam:
            if beam:
                # a target position has been skipped, so
                # add the cost of deleting the last character
                tdel_cost = self.tdel_cost(tpos-1) or 1.0e-20
                newcost -= numpy.log(tdel_cost)
            # make a new connection alignment
            new_alignment = self.copy(newcost=newcost)
            new_alignment.update(tpos)
            if verbosity:
                print(" New alg {}, cost {}".format(tuple(new_alignment), new_alignment.cost))
            alignments.append(new_alignment)
            tpos += 1
            beam += 1
        # make a new deletion alignment
        new_alignment = self.copy()
        new_alignment.update(-1)
        if verbosity:
            print(" New alg {}, cost {}".format(tuple(new_alignment), new_alignment.cost))
        alignments.append(new_alignment)
        for a in alignments:
            if a.complete():
                # complete alignment but we need to add the cost
                # of any remaining target chars
                last_tpos = a.last_nondel()
                tpos = last_tpos + 1
                newcost = 0
                additional = 0
                if verbosity:
                    print(" {} complete, last tpos {}".format(tuple(a), last_tpos))
                while tpos < self.tlength:
                    tdel_cost = self.tdel_cost(tpos) or 1.0e-20
                    log_td_cost = -numpy.log(tdel_cost)
                    if verbosity:
                        print("   tdel cost for {}: {}".format(tpos, log_td_cost))
                    newcost += log_td_cost
#                    new_alignment.cost -= numpy.log(tdel_cost)
                    tpos += 1
                    additional += 1
                a.cost += newcost
                if verbosity:
                    print("  {} extra tc, new cost {}, cost {}".format(additional, newcost, a.cost))
        return alignments

    ### Evaluation during search

    def prob1(self, sindex):
        """
        Probability associated with the connection in source
        word.
        """
        s = self.swi(sindex)
        # aligned position
        a = self[sindex]
        if a < 0:
            # character is deleted; get the del prob for s
            p = self.aligner.get_prob(s, -1)
        else:
            t = self.twi(a)
            p = self.aligner.get_prob(s, t)
        return p

    def prob(self):
        """
        Product of probabilities for current connections.
        """
        return product(*[self.prob1(i) for i in range(len(self))])

    def cost(self):
        """
        Negative log of probability product.
        """
        p = self.prob() or 1.0e-20
        return - numpy.log(p)

    def tdel_cost(self, tpos):
        """
        Get the cost the probability of
        deleting the target character at position tpos.
        """
        twi = self.twi(tpos)
        return self.aligner.tdeltas[twi]

    def distance(self):
        """
        Number of schars not yet connected.
        """
        return self.slength - len(self)

    def value(self):
        """
        Distance to end of word + negative log prob.
        """
        return self.distance() + self.cost()

    ### Evaluation of complete alignment, based on
    ### constraints

    def reverse(self):
        """
        The target-to-source alignment corresponding
        to this alignment.
        """
        tsa = []
        for sp, a in enumerate(self):
            # get the non-deleted alignment covered
            nondel = [p for p in self[:sp] if not self.deleted(p)]
            if nondel:
                last_tsa = nondel[-1]
            else:
                last_tsa = -1
            if a > last_tsa + 1:
                # one or more deleted characters in target
                tstart = last_tsa + 1
                for tp in range(tstart, a):
                    tsa.append(-1)
            if not self.deleted(a):
                # source char is aligned with something
                tsa.append(sp)
        # Positions in target beyond position aligned with end of source
        # last aligned position that is not -1
        if self.complete():
            last_nondel = self.last_nondel()
            if last_nondel < self.tlength:
                for ti in range(last_nondel+1, self.tlength):
                    tsa.append(-1)
        return tsa

    def recall(self, verbosity=0):
        """
        Return a score for this alignment based on
        constraints in the self.constraints list.
        """
        abs_scores = []
        int_scores = []
        if self.constraints:
            if not self.left:
                alignment = self.right2left()
            else:
                alignment = self
            # constraints is a pair:
            #  (absolute_constraints, interval_constraints)
            abs_c, int_c = self.constraints
            if abs_c:
                if verbosity:
                    print("Absolute constraints")
                score = 0
                total = 0
                for sp, tp in abs_c:
                    if verbosity:
                        print(" {}->{};".format(sp, tp), end='')
                    sat = False
                    sa = alignment[sp]
                    total += 1
                    if sa == tp:
                        sat = True
                        score += 1
                    if verbosity:
                        if sat:
                            print("  satisfied")
                        else:
                            print("  not satisfied")
                if verbosity:
                    print(" raw: {}".format(score))
                abs_scores.append((score, total))
#                    abs_score = score / total
            if int_c:
                a_rev = self.reverse()
                if verbosity:
                    print("Interval constraints")
                for beginning, end in int_c:
                    if verbosity:
                        print(" {}...{}->{}...{};".format(beginning[0], end[0], beginning[1], end[1]), end='')
                    good = 0
                    sp0, tp0 = beginning
                    sp1, tp1 = end
                    sln = sp1 - sp0 + 1
                    tln = tp1 - tp0 + 1
                    # score number of characters in longer segment
                    # which are inside the interval
                    if sln >= tln:
                        # check source characters
                        for sp in range(sp0, sp1+1):
                            ta = self[sp]
                            total += 1
                            if tp0 <= ta <= tp1:
                                good += 1
                            elif ta >= 0:
                                good -= 1
                        missed = tln - good
                        score = max([0, sln - missed])
                        int_scores.append((score, sln))
                    else:
                        # check target characters
                        for tp in range(tp0, tp1+1):
                            sa = a_rev[tp]
                            if sp0 <= sa <= sp1:
                                good += 1
                            elif sa >= 0:
                                good -= 1
                        missed = sln - good
                        score = max([0, tln - missed])
                        int_scores.append((score, tln))
                    if verbosity:
                        g, l = int_scores[-1]
                        print("  {} chars out of {}".format(g, l))
        return abs_scores, int_scores

    ### Displaying the alignment

    def pretty(self, verbosity=0, printit=False):
        """
        Pretty print alignment.
        """
        t_string = []
        s_string = []
        last_align = -1
        if self.left:
            schars = self.schars
            tchars = self.tchars
        else:
            schars = list(reversed(self.schars))
            tchars = list(reversed(self.tchars))
        for si, schar in enumerate(schars[:len(self)]):
            current_max = len(schar) + 1
            string = " " + schar
            align = self[si]
            if self.deleted(align):
                # deleted char in source
                t_string.append(' ' + ' ' * len(schar))
            elif align > last_align + 1:
                # one or more deleted characters in target
                tstart = last_align + 1
                for tindex in range(tstart, align):
                    tchar = tchars[tindex]
                    new_length = len(tchar) + 1
                    current_max = max(current_max, new_length)
                    s_string.append(" " + ' ' * len(tchar))
                    t_string.append(' ' + tchar)
            if not self.deleted(align):
                # source char is aligned with something
                last_align = align
                tchar = tchars[align]
                new_length = len(tchar) + 1
                current_max = max(current_max, new_length)
                newtstring = ''
                if si > 0 and align == self[si-1]:
                    # This character already part of string
                    newtstring = ' ' + ' ' * len(tchar)
                else:
                    newtstring = ' ' + tchar
                newtstring = newtstring.ljust(current_max)
                t_string.append(newtstring)
            string = string.ljust(current_max)
            s_string.append(string)
        # Positions in target beyond position aligned with end of source
        # last aligned position that is not -1
        if self.complete():
            last_nondel = self.last_nondel()
            if verbosity:
                print(" end of s, last t nondel {}".format(last_nondel))
            if last_nondel < self.tlength:
                for ti in range(last_nondel+1, self.tlength):
                    tchar = tchars[ti]
                    t_string.append(" " + tchar)
                    s_string.append(" " + " " * len(tchar))
        if not self.left:
            s_string.reverse()
            t_string.reverse()
        s_string = ''.join(s_string)
        t_string = ''.join(t_string)
        if printit:
            print(s_string)
            print(t_string)
        else:
            return s_string + "\n" + t_string
