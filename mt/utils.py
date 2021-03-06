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

-- 12.7.2020
   Created.
"""

import itertools, datetime

def get_file_id(filename):
    """
    Get the characters before the first number in the filename.
    """
    name = ''
    for char in filename:
        if char.isdigit():
            return name
        else:
            name += char
    return name

### Time formatting

TIME_FORMAT = "%d.%m.%Y/%H:%M:%S:%f"
# Time format without microseconds;
SHORT_TIME_FORMAT = "%d.%m.%Y/%H:%M:%S"
# Format without punctuation, for filenames (sortable by date created)
FILE_SHORT_TIME_FORMAT = "%Y%m%d%H%M%S"

ZERO_TIME = datetime.timedelta()
TIME0 = datetime.datetime.utcnow()

def get_time(string=False, short=False, file=False):
    time = datetime.datetime.utcnow()
    if string:
        return time2str(time, short=short, file=file)
    return time

def get_time_since0(time):
    return time - TIME0

def time2str(time, short=False, file=False):
    if file:
        format = FILE_SHORT_TIME_FORMAT
    elif short:
        format = SHORT_TIME_FORMAT
    else:
        format = TIME_FORMAT
    return time.strftime(format)

def str2time(string, short=False, file=False):
    if file:
        format = FILE_SHORT_TIME_FORMAT
    elif short:
        format = SHORT_TIME_FORMAT
    else:
        format = TIME_FORMAT
    return datetime.datetime.strptime(string, format)

### Combinations and probabilities (used in Aligner)

def product(*numbers):
    """
    Product of the arguments.
    """
    if not numbers:
        return 1
    else:
        return numbers[0] * product(*numbers[1:])

def seq_prob_indices(seq, indices):
    """
    Seq is a list of probabilities. Indices indicate
    a combination of elements from seq.
    Returns the product of the probabilities at indices
    and the complements of the probabilities in other positions.
    """
    probs = [(p if i in indices else 1-p) for i, p in enumerate(seq)]
    return product(*probs)

def seq_prob_comb(seq, n, icombs=None):
    """
    Seq is a list of probabilities. n is an integer specifying
    the number of probabilities to be sampled.
    Returns the sum of the products of the probabilities for
    every combination of n elements in seq.
    """
    combs = icombs or itertools.combinations(range(len(seq)), n)
    # print("  Combs for seq {}, n {}".format(seq, n))
    # print("  {}".format(list(combs)))
    products = [seq_prob_indices(seq, c) for c in combs]
    return sum(products)

def all_index_combs(n):
    """
    List of lists all combinations of integers in range 1 to n,
    one sublist for each combination length.
    """
    indices = range(n)
    return [list(itertools.combinations(indices, i+1)) for i in indices]

def weighted_seq_prob_combs(seq, icombs=None):
    """
    Seq is a list of probabilities associated with elements
    in another sequence of the same length.
    Returns the expected number of elements from the sequence.
    For every n in range(len(n)), use seq_prob_comb
    to calculate the sum of the products of the probabilities
    for every combination of n elements in seq.
    Then sum these sums, weighted by n.
    """
    if not icombs:
        icombs = all_index_combs(len(seq))
    weighted_combs = [(i+1) * seq_prob_comb(seq, i+1, icombs[i])\
                     for i in range(len(seq))]
    return sum(weighted_combs)
