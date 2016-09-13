#!/usr/bin/env python
"""
"""

import re
import ctypes
import warnings
from pkg_resources import resource_filename

from numpy.ctypeslib import load_library, ndpointer
import numpy as np

libRONN = None
"""C library that does the heavy-lifting (loaded upon 1st use)
"""

re_remove = None
"""Initialize regex used to remove invalid characters
"""

disorder_weight = 0.53
"""
"""

def load_libRONN(weight=disorder_weight):
    global libRONN
    libRONN = load_library("libRONN", resource_filename(__name__, '.'))

    # Read all data files for models
    libRONN.read_all_models.argtypes = [ctypes.c_char_p, ctypes.c_float]
    libRONN.read_all_models.restype = ctypes.c_int
    libRONN.read_all_models(resource_filename(__name__, './data'), weight)

    # Initialize predict_seq method
    libRONN.predict_seq.argtypes = [ctypes.c_char_p]
    libRONN.predict_seq.restype = ndpointer(dtype=np.float)

    return


def calc_ronn(seq):
    """Calculate hotloop, coil, and REM465 propensities

    Wraps predict_seq from libdisembl.c. Loads library upon first call.

    Parameters
    ----------
    seq : str
        Protein sequence for which propensities will be calculated

    handle_invalid : bool
        If `True`, attempt to handle invalid residues.

    Returns
    -------
    pd.DataFrame
        coils, hotloops, & rem465 propensities

    Raises
    ------
    None
    """

    if libRONN is None:
        load_libRONN()

        global re_remove
        re_remove = re.compile('[^FIVWMLCHYAGNRTPDEQSK]')

    # # seq_len = len(seq)
    # if len(seq) <= 20:
    #     # seq = "A" * 20 + seq + "A" * 20
    #     warnings.warn(("%s...%s: Sequence shorter than Neutral network window."
    #                   " Effectively padding with K") % (seq[:7], seq[-7:]))

    # # handle invalid characters
    # missing_pos = [match.span()[1]-1 for match in re_remove.finditer(seq)]
    # seq = re_remove.sub('', seq)


    # http://stackoverflow.com/a/37888716/2320823
    return libRONN.predict_seq(str.encode(seq))

    # if missing_pos and handle_invalid:
    #     warnings.warn("%s...%s: Unknown residues. Interpolating." %
    #                   (seq[:7], seq[-7:]))
    #     # If there are multiple missing, doing this in order will keep
    #     # indicies correct
    #     for pos in missing_pos:
    #         seq = seq[:pos] + 'X' + seq[pos:]
    #         try:
    #             # Average the values that are immediately surrounding
    #             coils = np.insert(coils, obj=pos,
    #                         values=np.mean((coils[pos-1], coils[pos])))
    #             rem465 = np.insert(rem465, obj=pos,
    #                         values=np.mean((rem465[pos-1], rem465[pos])))
    #             hotloops = np.insert(hotloops, obj=pos,
    #                         values=np.mean((hotloops[pos-1], hotloops[pos])))
    #         except IndexError as err:
    #             if pos == 0:
    #                 coils = np.insert(coils, obj=pos, values=coils[pos])
    #                 rem465 = np.insert(rem465, obj=pos, values=rem465[pos])
    #                 hotloops = np.insert(hotloops, obj=pos,
    #                                      values=hotloops[pos])
    #             elif pos == coils.size:
    #                 coils = np.append(coils, values=coils[-1])
    #                 rem465 = np.append(rem465, values=rem465[-1])
    #                 hotloops = np.append(hotloops, values=hotloops[-1])
    #             else:
    #                 print(pos, coils.size, flush=True)
    #                 raise IndexError(err)
    #
    # return pd.DataFrame({'residue': list(seq),
    #                      'coils': coils,
    #                      'rem465': rem465,
    #                      'hotloops': hotloops})
