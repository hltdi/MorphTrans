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

-- 29.1.2020
   Created.
"""

import os

from .rule import *
from .align import *

class Scoring:
    """
    Functions for figuring the cost of insertion, deletion,
    substitution, and no change in editing word.
    """

    CONTEXT_DECAY = [1.0, 0.8, 0.6, 0.25, 0.12]

    @staticmethod
    def context_cost(phone, context, forward=False):
        if phone not in context:
            return 0
        if forward:
            distance = context.index(phone)
        else:
            distance = len(context) - context.index(phone) - 1
        if distance > len(Scoring.CONTEXT_DECAY):
            return 0
        return Scoring.CONTEXT_DECAY[distance]

    @staticmethod
    def delete(basic_cost, length_cost, phone_cost):
        """Given basic deletion cost, cost for short source words,
        and the cost of deleting the particular phone,
        return a function that takes the phone to be deleted and
        difference in source and target word lengths."""
        def delhelp(phone, ldiff):
            cost = basic_cost
            if ldiff <= 0:
                cost += length_cost
            if phone_cost:
                cost += phone_cost * Phone.insdel_cost(phone)
            return cost
        return delhelp

    @staticmethod
    def insert(basic_cost, length_cost, context_cost, phone_cost,
               forward=False):
        """Given basic insertion cost, cost for long source words,
        cost for inserting a character that's already nearby in the
        source word, and the cost for inserting the particular phone,
        return a function that takes the phone to be inserted,
        the difference in source and target word lengths, and the
        (previous) source word context."""
        def inshelp(phone, ldiff, context, forward):
            cost = basic_cost
            if phone in context:
                # Penalize insertion of consonant that is in source context
                context_mult = Scoring.context_cost(phone, context, forward=forward)
                cost += context_mult * context_cost
            if ldiff >= 0:
                cost += length_cost
            if phone_cost:
                cost += phone_cost * Phone.insdel_cost(phone)
            return cost
        return inshelp

    @staticmethod
    def substitute(basic_cost, scontext_cost, tcontext_cost, phone_cost,
                   forward=False):
        """Given basic substitution cost, the cost for replacing a source
        phone that's already nearby in the source word, and the cost for
        replacing it with a target phone that nearby in the target word,
        return a function that takes the two phones and the two contexts."""
        def subshelp(newphone, oldphone, scontext, tcontext, forward):
            cost = basic_cost
            if newphone in scontext:
                context_mult = Scoring.context_cost(newphone, scontext, forward=forward)
#                print(">>newphone {} in scontext {}, mult {}".format(newphone, scontext, context_mult))
                cost += context_mult * scontext_cost
            if oldphone in tcontext:
                context_mult = Scoring.context_cost(oldphone, tcontext, forward=forward)
                cost += context_mult * tcontext_cost
            if phone_cost:
                distance = Phone.distance(newphone, oldphone)
#                print(">>distance between {} and {}: {}".format(newphone, oldphone, distance))
                cost += phone_cost * distance
            return cost
        return subshelp

    @staticmethod
    def nochange(basic_cost, match_cost):
        """Given the basic cost for making no change (substituting the
        source phone for itself) and the (negative) cost based on the number
        of phone matches in other target words in positions aligned with this
        source phone, return a function that takes the number of matches."""
        def nochhelp(matches):
            return basic_cost + match_cost * matches
        return nochhelp

class CRFProblem:
    """A dataset and scoring functions to use in learning rules for
    word-to-word translation based on the approach in Durrett & DeNero."""

    # default scoring functions
    default_scoring = {'insert': Scoring.insert(0, 0.25, 0.5, 1.0),
                       'delete': Scoring.delete(0, 0.25, 1.0),
                       'substitute': Scoring.substitute(0, 0.75, 0.75, 1.0),
                       'nochange': Scoring.nochange(0, -1.0)}

    id = 0

    def __init__(self, datafile='', data=None, info=None, scoring=None,
                 name=''):
        if data:
            # The data is passed to the constructor.
            self.data = data
            for tg in self.data:
                tg.problem = self
            self.info = info or {}
        else:
            # The data is read in from a file.
            datafile = datafile or '1.dt'
            self.data, self.info = Data.read(datafile, self)
            for i, tg in enumerate(self.data):
                self.data[i] = TGroup([Word(w) for w in tg['words']],
                                      problem=self)
        # A dict of key, functions for each
        # of delete, insert, substitute, and nochange.
        self.scoring = scoring or CRFProblem.default_scoring
        self.set_name(name)

    def set_name(self, name):
        if name:
            self.name = "Problem: {}".format(name)
        else:
            self.name = "Problem: {}".format(CRFProblem.id)
            CRFProblem.id += 1

    def __repr__(self):
        return self.name

    ## Alignment of TGroups in self.data

    def align(self, direction=None, verbosity=0):
        """Align each of the TGroups in the Problem. If direction is specified ('forward' or 'backward'),
        use this direction always. If not, try both directions for each TGroup, selecting the cheaper."""
        print("ALIGNING GROUPS in {}".format(self))
        total_cost = 0
        for i, tg in enumerate(self.data):
            print()
            if direction:
                print("ALIGNING TGROUP {}: {}, DIRECTION: {}".format(i, tg, direction))
                cost = tg.align(direction=direction, verbosity=verbosity)
            else:
                print("ALIGNING TGROUP {}: {}, BOTH DIRECTIONS".format(i, tg))
                # Try both directions for TGroup, selecting one with the lower cost
                forward_alignments = tg.copy_alignments()
                backward_alignments = tg.copy_alignments()
                forward_cost = tg.align(direction='forward', alignments=forward_alignments,
                                        verbosity=verbosity)
                backward_cost = tg.align(direction='backward', alignments=backward_alignments,
                                         verbosity=verbosity)
                if forward_cost < backward_cost:
                    cost = forward_cost
                    tg.alignments = forward_alignments
                    print("FORWARD MINIMIZATION CHEAPER FOR TGROUP {}".format(tg))
                else:
                    cost = backward_cost
                    tg.alignments = backward_alignments
                    print("BACKWARD MINIMIZATION CHEAPER FOR TGROUP {}".format(tg))
                print("FINAL ASSIGNMENTS:")
                tg.print_alignments()
            total_cost += cost
        print()
        print("TOTAL COST: {}".format(total_cost))

    def print_alignments(self, verbosity=0):
        for i, tg in enumerate(self.data):
            print("TGroup {}: {}".format(i, tg))
            tg.print_alignments(verbosity=verbosity)

    ## Functions for scoring
    def score_insert(self, index, phone, ldiff, context, forward=False):
        insert = self.scoring['insert']
        return insert(phone, ldiff, context, forward=forward)

    def score_delete(self, index, phone, ldiff):
        delete = self.scoring['delete']
        return delete(phone, ldiff)

    def score_compare(self, index, newphone, oldphone, scontext, tcontext, matches=0, forward=False):
        if newphone == oldphone:
            # More sophisticated comparison can happen here
#            print("SAME PHONE: {}, matches {}".format(newphone, matches))
            nochange = self.scoring['nochange']
            return nochange(matches)
        else:
            substitute = self.scoring['substitute']
            return substitute(newphone, oldphone, scontext, tcontext, forward=forward)
