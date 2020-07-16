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

from .data import *
from .word import *
from .utils import *
from .search import *
import numpy

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
    TRAINING = 0
    TEST = 1
    VALIDATION = 2

    def __init__(self, source, target, datafiles=None,
                 test=True, validate=True):
        self.source = source
        self.target = target
        # Source and target language characters
        self.schars = self.get_schars()
        self.tchars = self.get_tchars()
        # Number of source and target characters
        self.ns = len(self.schars)
        self.nt = len(self.tchars)
        # Character association probabilities
        self.probs = self.make_probs()
        self.init_probs()
        self.tnorms = numpy.zeros(self.nt+1)
        # History of sum of squares of changes in probs during M-step
        self.changes = []
        # Justification/alignment parameters
        self.lambdas = self.make_ljust_probs()
        ## Variables depending on data
        # Subsets of data for training, test, validation
        # List of word pairs, consisting of character strings
        self.data = []
        # List of word pairs, consisting of lists of character indices
        self.data_indices = []
        self.training = []
        self.training_indices = []
        self.test = []
        self.test_indices = []
        self.validation = []
        self.validation_indices = []
        # Other arrays that depend on data
        self.counts = []
        self.tchar_counts = []
        self.tdeltas = []
        self.tdelt_counts = []
        self.ljustcounts = []
        self.rjustcounts = []
        self.diffcounts = []
        self.tindex_combs = []
        if datafiles:
            for datafile in datafiles:
                d = Data(datafile)
                wpairs = d.get_words([self.source, self.target])
                self.data.extend(wpairs)
            # Shuffle data so that it can be split into training, testing
            # and validation subsets
            self.shuffle(self.data)
            if test or validate:
                self.test, self.validation, self.training =\
                self.split_data(self.data)
            else:
                self.training = self.data
            # Or only work with indexed data
            self.data_indices = [self.wpair2indices(wp) for wp in self.data]
            if test:
                self.test_indices = [self.wpair2indices(wp) for wp in self.test]
            if validate:
                self.validation_indices = [self.wpair2indices(wp) for wp in self.validation]
            self.training_indices = [self.wpair2indices(wp) for wp in self.training]
            self.trainingN = len(self.training)
            self.counts = self.make_counts()
            # Target character deletion parameters and counts
            self.tdelt_counts = self.make_tdel_counts()
            # Number of target characters
            self.tchar_counts = self.make_tchar_counts()
            self.tdeltas = self.make_tdeltas()
            # Store combinations of target indices to avoid
            # repeating itertools.combinations()
            self.tindex_combs = self.make_target_index_combs()
            # Alignment/justification counts
            self.ljustcounts = self.make_just_counts()
            self.rjustcounts = self.make_just_counts()
            self.diffcounts = self.make_diff_counts()

    def __repr__(self):
        return "↙↘:{},{}".format(self.source, self.target)

    ### datasets

    def shuffle(self, data, reproduce=True):
        """
        Shuffle items in dataset so they can be split into training, test,
        validation sets.
        """
        if reproduce:
            numpy.random.seed(Aligner.seedgen)
        numpy.random.shuffle(data)

    def split_data(self, data, test=True, validate=True):
        """Split data into test, validation, and training sets."""
        n = len(data)
        test_i = 0
        validation_i = 0
        if test:
            test_i = round(n * Aligner.test_frac)
        if validate:
            validation_i = test_i + round(n * Aligner.validation_frac)
        return data[:test_i], data[test_i:validation_i], data[validation_i:]

    def get_data_indices(self,
                         training=False, test=False, validate=False):
        """
        Word pairs as indices into various parameter arrays.
        """
        if training:
            return self.training_indices
        elif test:
            return self.test_indices
        elif validate:
            return self.validation_indices
        else:
            return self.data_indices

    def get_data_chars(self,
                       training=False, test=False, validate=False):
        """
        Word pairs as lists of characters.
        """
        if training:
            return self.training
        elif test:
            return self.test
        elif validate:
            return self.validation
        else:
            return self.data

    def get_data(self,
                 training=False, test=False, validate=False):
        """
        (data_indices, data_chars) for particular datasets
        """
        return \
        (self.get_data_indices(training=training, test=test, validate=validate),
         self.get_data_chars(training=training, test=test, validate=validate))

    ### Accessing characters, probabilities, indices

    def get_schars(self, every=False):
        phones = Phone.PHONES.get(self.source)
        if every:
            return phones[0] + phones[1]
        return phones[0]

    def get_tchars(self, every=False):
        phones = Phone.PHONES.get(self.target)
        if every:
            return phones[0] + phones[1]
        return phones[0]

    def get_tchar_index(self, char, check=False):
        if check:
            if char not in self.tchars:
                return -1
        return self.tchars.index(char)

    def get_schar_index(self, char, check=False):
        if check:
            if char not in self.schars:
                return -1
        return self.schars.index(char)

    def get_schar(self, index):
        if index >= self.ns:
            return Aligner.NULL
        return self.schars[index]

    def get_tchar(self, index):
        if index >= self.nt:
            return Aligner.NULL
        return self.tchars[index]

    def get_prob(self, sindex, tindex):
        return self.probs[sindex][tindex]

    def get_schar_probs(self, char):
        charindex = self.get_schar_index(char)
        return self.probs[charindex]

    def get_tchar_probs(self, char):
        charindex = self.get_tchar_index(char)
        return self.probs[:,charindex]

    def get_char_prob(self, schar, tchar):
        sindex = self.get_schar_index(schar)
        tindex = self.get_tchar_index(tchar)
        return self.get_prob(sindex, tindex)

    def get_target_probs(self, tword, sindex):
        """
        Given a sindex, find the st probabilities of all characters
        in tword, with NULL in last position.
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

    # def get_del_prob(self, sword, tword):
    #     """
    #     Get the probability that a character in source word will be
    #     deleted (connected to the NULL character).
    #     """
    #     m = len(sword)
    #     l = len(tword)
    #     count = self.deltas[m-l]
    #     return count / m

    def get_ljust_prob(self, sword, tword):
        m = len(sword)
        l = len(tword)
        index = self.get_lendiff_index(m-l)
        return self.lambdas[index]

    # def get_delta(self, lendiff):
    #     index = self.get_lendiff_index(lendiff)
    #     return self.deltas[index]

    def get_NULL_prob(self, si):
        return self.probs[si,self.nt]

    def get_lendiff(self, index):
        """Return the len(s) - len(t) difference associated with the
        count or parameter in the index position."""
        return Aligner.lendiff0 + index

    def get_lendiff_index(self, diff):
        """
        Return the index associated with the len(s) - len(t) difference
        in the count or parameter table.
        """
        return diff - Aligner.lendiff0

    # def get_wpair_lendiff(self, wpair):
    #     """
    #     Return the difference in lengths between the source and target
    #     words in wpair.
    #     """
    #     return len(wpair[0]) - len(wpair[1])

    def wpair2indices(self, wpair):
        """
        Convert words in wpair to lists of indices into arrays.
        """
        sindices = [self.get_schar_index(c) for c in wpair[0]]
        tindices = [self.get_tchar_index(c) for c in wpair[1]]
        return sindices, tindices

    def targetN(self):
        """
        Number of target probs and counts, with the NULL position
        at the end.
        """
        return self.nt+1

    ### Arrays and parameters

    def make_probs(self):
        """
        Make the 'st' table of probabilities: p(s|t) for s/t character
        combinations. source is rows, target columns.
        """
        array = numpy.full((self.ns, self.targetN()), Aligner.otherprob)
        return array

    def init_probs(self):
        """
        Initialize the probabilities in the table.
        """
        array = self.probs
        getrowchar = self.get_schar
        getcolindex = self.get_tchar_index
        for i, row in enumerate(array):
            char = getrowchar(i)
            colindex = getcolindex(char, True)
            if colindex >= 0:
                array[i,colindex] = Aligner.eqprob
            # Deletion probability
            array[i,self.nt] = 0.2
        # Normalize
        self.norm_cols()

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
        combinations. source is rows, target columns.
        """
        array = numpy.zeros((self.ns, self.targetN(), self.trainingN))
        return array

    def make_just_counts(self):
        """
        Make the left- or right-justified table of normalized counts, one for each
        combination of len(s) - len(t) distance."""
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

    # def make_deltas(self):
    #     """
    #     Make the parameters governing how many characters in a source
    #     sentence will be deleted (connected to the NULL character), given
    #     different differences between source and target word lengths.
    #     The numbers represent additional deleted characters beyond what
    #     must be the case because of the diffence in word lengths.
    #     """
    #     return numpy.array(
    #     [0.0, 0.0, 0.01, 0.05, 0.5, 0.0, 0.3, 0.5, 0.7, 1.0, 1.3,
    #      1.5, 2.0, 2.3, 2.5, 3.0, 4.3]
    #      )

    # def make_del_counts(self):
    #     """
    #     Initialize the array where expected counts of deleted characters
    #     are stored for each len(s) - len(t) difference.
    #     """
    #     array = numpy.zeros((Aligner.nlendiffs))
    #     return array

    def make_tdeltas(self):
        """
        Make the parameters governing whether characters in a
        target word will be deleted.
        """
        # array = numpy.full((self.nt), 0.1)
        array = numpy.zeros((self.nt))
        for tindex in range(self.nt):
            if self.tchar_counts[tindex] > 0:
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
        occurs in the training set.
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

    # def init_deltas(self):
    #     """
    #     Initialize the deltas (character deletion parameters).
    #     """
    #     self.deltas = self.make_deltas()

    ### EM

    def EM(self, iter_cutoff=10, error_cutoff=-10.0, verbosity=0):
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
        Do whatever initialization is required before the E-step on each
        EM iteration.
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
        #schar = self.get_schar(si)
        #tchar = self.get_tchar(ti)
        leftprox, rightprox = self.calc_cooccs(sword, si, tword, ti)
        if leftprox or rightprox:
            cooccs = (leftprox + rightprox) / 2.0
            self.counts[si, ti, wpindex] = (prob * cooccs) / denom
            # also normalize these?
            self.ljustcounts[diffindex] += (prob * leftprox)
            self.rjustcounts[diffindex] += (prob * rightprox)
        #return leftprox, rightprox

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
        # if schar not in sword or tchar not in tword:
        #     return 0.0, 0.0
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
#        return sword.count(schar) * tword.count(tchar)

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
        # if verbosity:
        #     print(" updating tdelts counts; tword {}, tdelts {}".format(tword, tdelts))
        for tindex, tdelt in zip(tword, tdelts):
            exptdel = tdelt * exptdels / tdeltsum
            self.tdelt_counts[tindex] += exptdel
            if verbosity > 1:
                print("  Updating tdcount for {}: {}".format(tindex,
                self.tdelt_counts[tindex]))

    # def get_deletions(self, diff):
    #     """
    #     Get the expected number of deleted characters in source word of
    #     length slength associated with target word of length tlength.
    #     """
    #     param = self.get_delta(diff)
    #     baseline = max([0, diff])
    #     return baseline + param
    #     #return (baseline + param) / slength
    #
    #     """
    #     Calculate the exptected number of source word
    #     character deletions, based on the target word
    #     character probabilities and the difference in
    #     word lengths.
    #     """
    #     ldiff = len(sword) - len(tword)
    #     tdels = self.calc_exp_tdels(tword)
    #     return max(0, ldiff + tdels)

    def calc_exp_tdels(self, tword, tdelts):
        """
        Calculate the expected number of target char
        deletions, based on their character deletion
        probabilities.
        """
        # if tdelts == None:
        #     tdelts = self.get_tdeltas(tword)
        icombs = self.tindex_combs[len(tword)]
        return weighted_seq_prob_combs(tdelts, icombs)

    def calc_exp_sdels(self, sword, tword, tdels):
        """
        Calculate the expected number of source char
        deletions, based on the expected target char deletions
        and the lengths of the words.
        """
        ldiff = len(sword) - len(tword)
        # if not tdels:
        #     tdels = self.calc_exp_tdels(tword)
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

    #     for lendiffi in range(Aligner.nlendiffs):
    #         # length difference for this parameter
    #         lendiff = self.get_lendiff(lendiffi)
    #         # required deletions
    #         baseline = max([0, lendiff])
    #         # total estimated deletions for this length difference
    #         count = self.delcounts[lendiffi]
    #         # number of word pairs with this length difference
    #         nwpairs = self.diffcounts[lendiff]
    #         # estimated deletions per word pair
    #         deletions = count / nwpairs
    #         # number of extra deletions (greater than length difference)
    #         deletions = deletions - baseline

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

    ### Alignments for word pairs.

    def make_alignment(self, data_i, explicit=None,
                       training=False, test=False, validate=False):
        """
        Create an Alignment object for the wpair at index data_i.
        """
        indices, data = self.get_data(training=training,
                                      test=test,
                                      validate=validate)
        return Alignment(indices[data_i], data[data_i],
                         self, explicit=explicit)

    def make_align_searcher(self, data_i, dataset=-1):
        """
        Create a best-first searcher for alignments.
        """
        training = test = validate = False
        if dataset == Aligner.TRAINING:
            training = True
        elif dataset == Aligner.TEST:
            test = True
        elif dataset == Aligner.VALIDATION:
            validate = True

        searcher = \
        BestFirst("AlignSearcher",
                  goal_test=lambda s: s.complete(),
                  make_start=lambda: self.make_alignment(data_i,
                       training=training, test=test, validate=validate),
                  extend=lambda s: s.extend(),
                  evaluate=lambda s: s.cost + s.distance)
        return searcher

    @staticmethod
    def argmaxes(array):
        """
        Utility function. Like argmax except that it returns a list of
        positions in case there's a tie.
        """
        indices = [0]
        maximum = array[0]
        for i, a in enumerate(array[1:]):
            if a > maximum:
                maximum = a
                indices = [i+1]
            elif a == maximum:
                indices.append(i+1)
        return indices

#     def find_best_alignment(self, wpair):
#         sword, tword = wpair
#         alignment = []
#         conflicts = []
#         for scindex, sindex in enumerate(sword):
#             tprobs = self.get_target_probs(tword, sindex)
#             maxtcindices = Aligner.argmaxes(tprobs)
#             if len(maxtcindices) > 1:
#                 # the best target char appears more than once;
#                 # pick the closest one
#                 print("Multiple maxtcindices for {}|{}: {}".format(scindex, sindex, maxtcindices))
#                 stdists = [(abs(scindex - tcindex), tcindex) for tcindex in maxtcindices]
#                 stdists.sort(key=lambda x: x[0])
#                 besttcindex = stdists[0][1]
#             else:
#                 besttcindex = maxtcindices[0]
#             print("Best tc index for {}|{}: {}".format(scindex, sindex, besttcindex))
#             if besttcindex in alignment:
#                 prevscindex = alignment.index(besttcindex)
#                 conflicts.append(((prevscindex, scindex), besttcindex))
# #                if self.keep_old_stconn(sword, tword, scindex, prevscindex, besttcindex):
#                     # Previous connection is better; connect this scindex to NONE
# #                    alignment.append(len(tword))
# #                else:
#                     # New connection is better; replace old one with NONE
# #                    alignment[prevscindex] = len(tword)
#                     #alignment.append(besttcindex)
#             alignment.append(besttcindex)
#         return alignment, conflicts

    # def keep_old_stconn(self, sword, tword, newsci, oldsci, tci):
    #     """
    #     When two source char indices are to be connected to the same
    #     target char index (other than NONE), pick the scindex that is closer
    #     in relative index to the target char index.
    #     """
    #     olddist = Aligner.get_rel_index_distance(sword, oldsci, tword, tci)
    #     newdist = Aligner.get_rel_index_distance(sword, newsci, tword, tci)
    #     return olddist < newdist

class Alignment(list):
    """
    Alignment between a source word and its translation in
    one target language.
    """

    beam = 2
    '''Number of succeeding positions to align with current
    s position when extending.'''

    def __init__(self, wpair, wpair_chars, aligner,
                 cost=0, distance=0, explicit=None):
        if explicit:
            list.__init__(self, explicit)
        else:
            list.__init__(self)
        self.sword = wpair[0]
        self.tword = wpair[1]
        self.slength = len(self.sword)
        self.tlength = len(self.tword)
        # We need this because it's got all of the parameters
        # and characters
        self.aligner = aligner
        self.schars = wpair_chars[0]
        self.tchars = wpair_chars[1]
        # Accumulated evaluation
        self.cost = cost
        self.distance = distance or len(self.sword)

    def __repr__(self):
        #return "{}→{}\n{}".format(
        #''.join(self.schars),
        #''.join(self.tchars),
        return self.pretty()

    ### Access

    def swi(self, index):
        return self.sword[index]

    def twi(self, index):
        return self.tword[index]

    def schar(self, index):
        return self.schars[index]

    def tchar(self, index):
        return self.tchars[index]

    def last_nondel(self):
        nondel = [p for p in self if not self.deleted(p)]
        if nondel:
            return nondel[-1]
        else:
            return -1

    ### Extending: creating new search states

    def copy(self):
        """
        Return a copy of the Alignment with the same
        word pair and indices.
        """
        return Alignment((self.sword, self.tword),
                         (self.schars, self.tchars),
                         self.aligner,
                         explicit=self[:],
                         cost=self.cost,
                         distance=self.distance)

    def update(self, tpos):
        """
        Update an alignment copy with the latest index.
        """
        self.append(tpos)
        self.update_cost()
        self.distance -= 1

    def deleted(self, position):
        return position < 0

    def complete(self):
        return len(self) == self.slength

    # def add_delete(self):
    #     self.append(-1)

    def next(self):
        """
        The character following the last one aligned.
        """
        if self.complete():
            # No more characters
            return
        return self.sword[len(self)]

    def extend(self):
        if self.complete():
            return
        alignments = []
        # make new alignments of the current source position
        # with the next beam positions
        if len(self) == 0 or all([self.deleted(p) for p in self]):
            last_tpos = -1
        else:
            last_tpos = self.last_nondel()
#            last_tpos = self[-1]
        tpos = last_tpos + 1
        beam = 0
        while tpos < self.tlength and beam < Alignment.beam:
            new_alignment = self.copy()
            new_alignment.update(tpos)
            alignments.append(new_alignment)
            tpos += 1
            beam += 1
        # make a new deletion alignment
        new_alignment = self.copy()
        new_alignment.update(-1)
        alignments.append(new_alignment)
        return alignments

    ### Evaluation

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

    def update_cost(self):
        """
        Add cost of last index.
        """
        index = len(self) - 1
        self.cost -= numpy.log(self.prob1(index))

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

    def pretty(self, verbosity=0, print=False):
        """
        Pretty print alignment.
        """
        t_string = []
        s_string = []
        last_align = 0
        for si, schar in enumerate(self.schars[:len(self)]):
            current_max = len(schar) + 1
            string = " " + schar
            align = self[si]
            if self.deleted(align):
                # deleted char in source
                t_string.append('  ')
            elif align > last_align + 1:
                # one or more deleted characters in target
                tstart = last_align + 1
                for tindex in range(tstart, align):
                    tchar = self.tchars[tindex]
                    new_length = len(tchar) + 1
                    current_max = max(current_max, new_length)
                    s_string.append(" " + ' ' * len(tchar))
                    t_string.append(' ' + tchar)
            if not self.deleted(align):
                # source char is aligned with something
                last_align = align
                tchar = self.tchars[align]
                new_length = len(tchar) + 1
                current_max = max(current_max, new_length)
                if si > 0 and align == self[si-1]:
                    # This character already part of string
                    t_string.append(' ' + ' ' * len(tchar))
                else:
                    t_string.append(' ' + tchar)
            string = string.rjust(current_max)
            s_string.append(string)
        # Positions in target beyond position aligned with end of source
        # last aligned position that is not -1
        if self.complete():
            last_nondel = self.last_nondel()
            if last_nondel < self.tlength:
                for ti in range(last_nondel+1, self.tlength):
                    t_string.append(" " + self.tchars[ti])
        s_string = ''.join(s_string)
        t_string = ''.join(t_string)
        if print:
            print(s_string)
            print(t_string)
        else:
            return s_string + "\n" + t_string
