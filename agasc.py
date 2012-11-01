#!/usr/bin/env python
from Chandra.Time import DateTime
import Ska.Shell
import Ska.Numpy


def agasc(ra, dec, radius=1.5, date=DateTime(), pm_correct=True):
    agasc_start_date = DateTime('2000:001:00:00:00.000')
    dsecs = date.secs - agasc_start_date.secs
    dyear = dsecs / (86400 * 365.25)
    milliarcsecs_per_degree = 3600 * 1000
    agasc_lines, denv = Ska.Shell.bash_shell(
        "/proj/sot/ska/bin/mp_get_agasc.pl -ra %s -dec %s -radius %s"
        % (ra, dec, radius))
    table = Ska.Table.read_ascii_table(agasc_lines)
    if not pm_correct:
        return table

    ra_corr = table.RA.copy()
    has_ra_pm = table['PM_RA'] != -9999
    ra_corr[has_ra_pm] = (table[has_ra_pm].RA
                          + (table[has_ra_pm].PM_RA
                             * (dyear / milliarcsecs_per_degree)))

    dec_corr = table.DEC.copy()
    has_dec_pm = table['PM_DEC'] != -9999
    dec_corr[has_dec_pm] = (table[has_dec_pm].DEC
                            + (table[has_dec_pm].PM_DEC
                                * (dyear / milliarcsecs_per_degree)))

    add_ra = Ska.Numpy.add_column(table, 'RA_PMCORR', ra_corr)
    corr = Ska.Numpy.add_column(add_ra, 'DEC_PMCORR', dec_corr)

    return corr
