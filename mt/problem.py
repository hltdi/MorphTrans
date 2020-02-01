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
from .data import *

class Problem:
    """A dataset and scoring functions to use in learning rules for
    word-to-word translation based on the approach in Durrett & DeNero."""

    # default scoring functions
    default_scoring = {'insert': Scoring.insert(1, 0.25, 0.5),
                       'delete': Scoring.delete(1, 0.25),
                       'substitute': Scoring.substitute(1, 0.5, 0.5),
                       'nochange': Scoring.nochange(0, -1.0)}

    def __init__(self, datafile='', data=None, info=None, scoring=None):
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
        # A dict of key, functions for each
        # of delete, insert, substitute, and nochange.
        self.scoring = scoring or Problem.default_scoring

    ## Functions for scoring
    def score_insert(self, index, phone, ldiff, context):
        insert = self.scoring['insert']
        return insert(phone, ldiff, context)

    def score_delete(self, index, phone, ldiff):
        delete = self.scoring['delete']
        return delete(phone, ldiff)

    def score_compare(self, index, newphone, oldphone, scontext, tcontext, matches=0):
        if newphone == oldphone:
            # More sophisticated comparison can happen here
            nochange = self.scoring['nochange']
            return nochange(matches)
        else:
            substitute = self.scoring['substitute']
            return substitute(newphone, oldphone, scontext, tcontext)

class Scoring:
    """
    Functions for figuring the cost of insertion, deletion,
    substitution, and no change in editing word.
    """

    @staticmethod
    def delete(basic_cost, length_cost):
        """Given basic deletion cost and cost for short source words,
        return a function that takes the phone to be deleted and
        difference in source and target word lengths."""
        def delhelp(phone, ldiff):
            if ldiff <= 0:
                lcost = length_cost
            else:
                lcost = 0
            return basic_cost + lcost
        return delhelp

    @staticmethod
    def insert(basic_cost, length_cost, context_cost):
        """Given basic insertion cost, cost for long source words,
        and cost for inserting a character that's already nearby in the
        source word, return a function that takes the phone to be inserted,
        the difference in source and target word lengths, and the
        (previous) source word context."""
        def inshelp(phone, ldiff, context):
            cost = basic_cost
            if phone in context:
                cost += context_cost
            if ldiff >= 0:
                cost += length_cost
            return cost
        return inshelp

    @staticmethod
    def substitute(basic_cost, scontext_cost, tcontext_cost):
        """Given basic substitution cost, the cost for replacing a source
        phone that's already nearby in the source word, and the cost for
        replacing it with a target phone that nearby in the target word,
        return a function that takes the two phones and the two contexts."""
        def subshelp(newphone, oldphone, scontext, tcontext):
            cost = basic_cost
            if newphone in scontext:
                cost += scontext_cost
            if oldphone in tcontext:
                cost += tcontext_cost
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

