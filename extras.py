import numpy as np
from astropy import log

# TODO: PR this over to pyspeckit maybe?
#       with complex guesses becoming available it's important to
#       pre-validate the guesses against parameter limits...
def verify_parlims(guesses, fiteach_args, npeaks, npars):
    """ Ensures that no guesses are over the min / max limits """
    runaways = []
    for idx in range(npeaks*npars):
        g_min, g_max = np.nanmin(guesses[idx]), np.nanmax(guesses[idx])
        minoverflow = g_min < fiteach_args['minpars'][idx]
        maxoverflow = g_max > fiteach_args['maxpars'][idx]
        zero_range = (fiteach_args['minpars'][idx]
                      == fiteach_args['maxpars'][idx])
        if zero_range and not fiteach_args['fixed'][idx]:
            log.error("Parameter #{} lower limit ({:.2f}) equals its higher"
                      " limit ({:.2f}). Set fixed[{}] to True?".format(
                      idx, fiteach_args['minpars'][idx],
                      fiteach_args['minpars'][idx], idx))
            runaways.append(idx)
        elif minoverflow:
            log.error("Parameter #{} minimum guess ({:.2f}) below"
                      " threshold ({:.2f})".format(idx,
                      g_min, fiteach_args['minpars'][idx]))
            runaways.append(idx)
        elif maxoverflow:
            log.error("Parameter #{} maximum guess ({:.2f}) above"
                      " threshold ({:.2f})".format(idx,
                      g_max, fiteach_args['maxpars'][idx]))
            runaways.append(idx)
    if runaways:
        raise RuntimeError("Stop! Fiteach will crash - parameter(s)"
                           " #{} out of bounds!".format(runaways))
