#!/usr/bin/env python
# Script to update guide star statistics


import os
import logging
import logging.handlers
import time
import mx.DateTime
import smtplib
from email.mime.text import MIMEText

import numpy as np
import Ska.DBI
import Ska.Table
import Ska.Numpy
import Ska.engarchive.fetch as fetch
import Ska.quatutil
import Ska.astro
from scipy.stats import scoreatpercentile
from Chandra.Time import DateTime
import agasc
from pexpect import ExceptionPexpect


sqlaca = Ska.DBI.DBI(dbi='sybase', server='sybase', user='aca_read',
                     numpy=True, database='aca')
sqlocc = Ska.DBI.DBI(dbi='sybase', server='sqlocc', user='aca_ops',
                     numpy=True, database='axafocat')


logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

ID_DIST_LIMIT = 1.5
acq_anom_radius = 160
mp_path = '/data/mpcrit1/mplogs'
revision = '3.3'

data_table = {'gui': 'trak_stats_data',
              'acq': 'acq_stats_data'}
warning_table = {'gui': 'trak_slot_warnings',
                 'acq': 'acq_stats_warnings'}

anom_list = os.path.join(os.environ['SKA'],
                         'data', 'star_stat_db', 'acq_anom_list.csv')
#anom_list = 'acq_anom_list.csv'


acq_cols = ['obsid',
            'obi',
            'tstart',
            'tstop',
            'slot',
            'cat_pos',
            'idx',
            'type',
            'agasc_id',
            'obc_id',
            'yang',
            'zang',
            'mag',
            'color',
            'halfw',
            'mag_obs',
            'yang_obs',
            'zang_obs',
            'd_mag',
            'd_yang',
            'd_zang',
            'y_offset',
            'z_offset',
            'ap_date',
            'revision']

gui_cols = ['obsid',
            'obi',
            'slot',
            'idx',
            'cat_pos',
            'id',
            'type',
            'color',
            'cyan_exp',
            'czan_exp',
            'mag_exp',
            'kalman_datestart',
            'kalman_datestop',
            'kalman_tstart',
            'kalman_tstop',
            'aoacyan_min',
            'aoacyan_mean',
            'aoacyan_max',
            'aoacyan_rms',
            'aoacyan_median',
            'aoacyan_5th',
            'aoacyan_95th',
            'aoaczan_min',
            'aoaczan_mean',
            'aoaczan_max',
            'aoaczan_rms',
            'aoaczan_median',
            'aoaczan_5th',
            'aoaczan_95th',
            'aoacmag_min',
            'aoacmag_mean',
            'aoacmag_max',
            'aoacmag_rms',
            'aoacmag_median',
            'aoacmag_5th',
            'aoacmag_95th',
            'n_samples',
            'not_tracking_samples',
            'bad_status_samples',
            'obc_bad_status_samples',
            'common_col_samples',
            'sat_pix_samples',
            'def_pix_samples',
            'quad_bound_samples',
            'ion_rad_samples',
            'mult_star_samples',
            'sample_interval_secs',
            'ap_date',
            'revision']

star_columns = [x for x in (set(acq_cols) | set(gui_cols))]


class TooNewError(Exception):
    """
    Special error for the case when the obsid is perhaps too new to have
    data in the starcheck or observations tables
    """
    pass


class TooOldError(Exception):
    """
    For observations before the observations table min tstart
    """
    pass


class WeirdObsidError(Exception):
    """
    For observations without star catalogs
    """
    pass


class MissingDataError(Exception):
    """
    Special error for the case when the obsid is missing from the
    starchec or observation database tables
    """
    pass


class ObsidError(Exception):
    """Class for obsid-related errors in this module"""
    pass


def get_options():
    from optparse import OptionParser
    parser = OptionParser(usage='update_guide_stats_db.py [options]')
    parser.set_defaults()
    parser.add_option("--dbi",
                      default='sqlite',
                      help="Database interface (sqlite|sybase)")
    parser.add_option("--server",
                      default='db_base.db3',
                      help="DBI server (<filename>|sybase)")
    parser.add_option("--user",
                      help="DBI user (Ska.DBI default)")
    parser.add_option("--database",
                      help="DBI database (Ska.DBI default)")
    parser.add_option("--obsid",
                      help="obsid to process \"manually\"")
    parser.add_option("--dryrun",
                      action='store_true',
                      )
    parser.add_option("--update-missing",
                      action="store_true",
                      help="add nominal no-starcat obsids to skip list")
    parser.add_option("--missing-list",
                      default="ok_missing.csv")
    parser.add_option('--email',
                      default="jeanconn",
                      help="email warning recipient")
    parser.add_option("-v", "--verbose",
                      type='int',
                      default=1,
                      help="Verbosity (0=quiet, 1=normal, 2=debug)")
    (opt, args) = parser.parse_args()
    return opt, args


def get_obs_db(obsid, obi, tstart):
    """
    Retrieve the observation table entry for the obsid/obi

    :param obsid: obsid (integer)
    :param obi: obi (integer)
    :rtype: recarray of database entry
    """

    obs_query = """select obsid, obi_num as obi,
                   kalman_datestart, kalman_datestop,
                   kalman_tstart, kalman_tstop
                   from observations
                   where obsid = %d and obi = %d""" % (obsid, obi)
    obs_array = sqlaca.fetchall(obs_query)
    if not len(obs_array):
        in_obs_all = sqlaca.fetchall("""select * from observations_all
                                        where obsid = %d and obi = %d"""
                                     % (obsid, obi))
        if len(in_obs_all):
            # this obsid doesn't have a kalman interval
            raise WeirdObsidError(
                "In observations_all but not observations table; no kalman")
        # bigger warning if this is old data and not in observations table
        minus7 = (DateTime(DateTime().mxDateTime
                           - mx.DateTime.DateTimeDeltaFromDays(7)))
        if (tstart >= minus7.secs):
            raise TooNewError("Not yet in observations table")
        else:
            # min kalman_tstart in observations table is 63073857.260169
            if (tstart < 63073857.260169):
                raise TooOldError("Missing from observations table; old obsid")
            raise MissingDataError("Missing from observations table")
    obsid = obs_array[0]
    return obsid


def get_needed_obsids(requested_obsid=None, missing_set=set()):
    """
    Fetch the list of to-be-done obsids.  If requested_obsid is
    specified, the details from that obsid will be the only entries in
    the returned structure.

    :param requested_obsid: optional obsid
    :rtype: list of recarrays with data from the mp_load_info table on obsid(s)
    """

    if requested_obsid:
        requested_obsid = int(requested_obsid)

    fields = ['mp.obsid', 'mp.obi',
              'mp.tstart', 'mp.tstop', 'mp.mp_path',
              'mp.last_ap_date', 'mp.no_starcheck', 'mp.wrong_starcheck']
    if (requested_obsid):
        query = """SELECT %s from mp_load_info as mp where obsid = %d""" % (
            ','.join(fields), requested_obsid)
        obsdata = sqlaca.fetchall(query)
        if not len(obsdata):
            logger.warn(
                "get_needed_obsids(): no record of %d in mp_load_info table"
                % requested_obsid)
    else:
        acq_q = ("""select %s, a.slot from mp_load_info as mp
               left outer join %s as a
               on mp.obsid = a.obsid
               and mp.obi = a.obi
               where a.slot is NULL
               order by mp.tstart desc"""
                 % (','.join(fields), data_table['acq']))
        #where a.slot = 4 or a.slot is NULL""" % (','.join(fields))
        # using the slot seems to be the quickest way to just get me one
        # entry per obsid...
        acq_list = sqlaca.fetchall(acq_q)
        acq_set = set((x['obsid'], x['obi']) for x in acq_list)

        gui_q = ("""select %s, a.slot from mp_load_info as mp
               left outer join %s as a
               on mp.obsid = a.obsid
               and mp.obi = a.obi
               where a.slot is NULL
               order by mp.tstart desc"""
                 % (','.join(fields), data_table['gui']))
        #where a.slot = 4 or a.slot is NULL""" % (','.join(fields))
        # using the slot seems to be the quickest way to just get me one
        # entry per obsid...
        gui_list = sqlaca.fetchall(gui_q)
        gui_set = set((x['obsid'], x['obi']) for x in gui_list)

        acq_up = acq_set - missing_set
        gui_up = gui_set - missing_set

        obsdata = []
        for t_obsid, t_obi in acq_up:
            obsdata.append(acq_list[(acq_list['obsid'] == t_obsid)
                                     & (acq_list['obi'] == t_obi)][0])

        for t_obsid, t_obi in gui_up - acq_up:
            obsdata.append(gui_list[(gui_list['obsid'] == t_obsid)
                                     & (gui_list['obi'] == t_obi)][0])

    return obsdata


def search_agasc(yang, zang, field_agasc, q_aca):
    """
    Search the retrieved agasc region for a star at the specified
    yang, zang, and return the star if there is a match.

    :param yang:
    :param zang:
    :param field_agasc: the retrieved agasc star field
    :param q_aca: pointing quaternion for obsid
    :rtype: recarray of the matching star or None
    """

    for agasc_star in field_agasc:
        ra, dec = Ska.quatutil.yagzag2radec(
            yang * 1.0 / 3600,
            zang * 1.0 / 3600,
            q_aca)
        # 3600*(sph_dist in degrees) for arcseconds
        dist = 3600 * Ska.astro.sph_dist(agasc_star['RA_PMCORR'],
                                         agasc_star['DEC_PMCORR'],
                                         ra, dec)
        if dist <= ID_DIST_LIMIT:
            return agasc_star

    return None


def get_stars(obsdb_obs, mp_obs, dbi):
    """
    Retrieve guide star catalog details from starcheck_catalog database and
    perform some basic checking (comparing stars to agasc and to archived
    stars). Return a dictionary of dictionaries, with slots as the top level
    keys. Each slot dictionary has a 'star' key and a 'warnings' key.  'star'
    corresponds to a recarray that is intended to be the database entry for
    the star. 'warnings' corresponds to a list of warnings for the slot.

    :param obs: basic obsid data as returned by get_needed_obsids()
    :rtype: dict
    """

    starcat = sqlaca.fetchall("""SELECT obsid,obi,slot,idx,id as agasc_id,
                                 idnote,type,yang,zang,mag,halfw as halfwidth
                                 from starcheck_catalog
                                 where obsid = %d and obi = %d
                                 order by idx"""
                              % (obsdb_obs.obsid, obsdb_obs.obi))

    soe_stars = sqlocc.fetchall("""SELECT slot,id,type,y_ang,z_ang from stars
                                   where obsid = %d
                                   and obi = %d order by slot"""
                                % (obsdb_obs.obsid, obsdb_obs.obi))

    if (len(starcat) and not len(soe_stars)):
        logger.warn("No SOE stars in archive for obsid %d obi %d"
                    % (obsdb_obs.obsid, obsdb_obs.obi))

    if not len(starcat):
        has_starcheck = sqlaca.fetchone("""SELECT obsid,obi from starcheck_obs
                                           where obsid = %d and obi = %d""" %
                                        (obsdb_obs.obsid, obsdb_obs.obi))
        if has_starcheck:
            raise WeirdObsidError(
                "In starcheck_obs but not starcheck_catalog; no catalog")
        if obsdb_obs.kalman_datestart < '2002:000:00:00:00.000':
            raise TooOldError("No starcheck_obs entry; old obsid")
        if obsdb_obs.obsid > 60000:
            raise WeirdObsidError("No catalog.  Expected weird obsid > 60000")

        raise ObsidError(
            "No guide stars found in starcheck_catalog for obsid %d obi %d"
            % (obsdb_obs.obsid, obsdb_obs.obi))

    starcheck_warnings = sqlaca.fetchall("""SELECT *
                                         from starcheck_warnings
                                         where obsid = %d and obi = %d """
                                         % (obsdb_obs.obsid, obsdb_obs.obi))

    obs_info = sqlaca.fetchone("""SELECT * from starcheck_obs
                                  where obsid = %d and obi = %d"""
                               % (obsdb_obs.obsid, obsdb_obs.obi))

    stars = np.rec.fromrecords([[None] * len(star_columns)] * len(starcat),
                               names=(star_columns))

    warnings = {}

    if not obs_info:
        return stars

    from Quaternion import Quat
    q_aca = Quat((obs_info['point_ra'],
                  obs_info['point_dec'],
                  obs_info['point_roll']))
    field_agasc = agasc.get_agasc_cone(ra=obs_info['point_ra'],
                                       dec=obs_info['point_dec'],
                                       radius=1.5,
                                       date=DateTime(obsdb_obs['kalman_tstart']))

    # position in catalog
    acq_cat_pos = 0
    for s in starcat:
        #if (stat_type == 'acq'):
        #    if ((s['type'] != 'ACQ') and (s['type'] != 'BOT')):
        #        continue

        # populate a recarray for the star, leaving the items expected from
        # telemetry empty
        star = stars[s['idx'] - 1]
        star.obsid = obsdb_obs.obsid
        star.obi = obsdb_obs.obi
        star.slot = s.slot
        star.idx = s.idx
        star.type = s.type
        if dbi == 'sqlite':
            star.ap_date = str(mp_obs.last_ap_date)
        else:
            star.ap_date = mp_obs.last_ap_date
        star.revision = "%s" % revision
        star.tstart = mp_obs.tstart
        star.tstop = mp_obs.tstop
        star.halfw = s.halfwidth
        star.yang = s.yang
        star.zang = s.zang
        star.mag = s.mag
        star.agasc_id = s.agasc_id
        star.kalman_datestart = obsdb_obs.kalman_datestart
        star.kalman_datestop = obsdb_obs.kalman_datestop
        star.kalman_tstart = obsdb_obs.kalman_tstart
        star.kalman_tstop = obsdb_obs.kalman_tstop
        star.cyan_exp = s.yang
        star.czan_exp = s.zang
        star.mag_exp = s.mag
        star.id = s.agasc_id
        star_warnings = []
        # reduced operations if this is a FID
        # (the storage operations are duplicated in the code)
        if s.type == 'FID':
            warnings[star.idx] = star_warnings
            continue

        if s['agasc_id'] is None:
            star_warnings.append("Missing AGASC id in starcheck_catalog")

        # for stars table 'type' field filtering
        soe_type = 1
        # increment the acquisition star catalog position
        # for the silly ACQID by catalog position stuff...
        if ((s['type'] == 'ACQ') or (s['type'] == 'BOT')):
            star.cat_pos = acq_cat_pos
            soe_type = 0
            acq_cat_pos += 1

        if len(soe_stars) == 0:
            star_warnings.append("SOE stars table has no entry for slot")
        else:
            soe_match = soe_stars[(soe_stars.slot == star.slot)
                                  & (soe_stars.type == soe_type)]
            if len(soe_match):
                # if there is a official entry for the star and it has
                # an agasc id, do they match?
                if s['agasc_id'] is not None:
                    if (soe_match['id'][0] != s['agasc_id']):
                        star_warnings.append("SOE id %d != Starcheck id %d" %
                                             (soe_match['id'][0],
                                              s['agasc_id']))
            else:
                star_warnings.append("SOE stars table has no entry for slot")

        # AGASC search
        agasc_match = 0
        # if the star that is in the starcheck catalog exists in the field
        if s['agasc_id'] and any(field_agasc['AGASC_ID'] == s['agasc_id']):
            agasc_star = field_agasc[
                field_agasc['AGASC_ID'] == s['agasc_id']][0]
            # does it have the same position as the starcheck entry?
            star_ra, star_dec = Ska.quatutil.yagzag2radec(
                s['yang'] * 1.0 / 3600,
                s['zang'] * 1.0 / 3600,
                q_aca)
            radial_dist = 3600 * Ska.astro.sph_dist(star_ra, star_dec,
                                                    agasc_star['RA_PMCORR'],
                                                    agasc_star['DEC_PMCORR'])
            if radial_dist < ID_DIST_LIMIT:
                agasc_match = 1
                star.mag = agasc_star['MAG_ACA']
                star.mag_exp = agasc_star['MAG_ACA']
                star.color = agasc_star['COLOR1']
            else:
                star_warnings.append(
                    "Starcheck ID not consistent with AGASC"
                    " position (off by %s arcsec)"
                    % radial_dist)

        # if we have no entry or starcheck doesn't match expectations
        if not agasc_match:
            agasc_star = search_agasc(star.cyan_exp,
                                      star.czan_exp,
                                      field_agasc,
                                      q_aca)
            if agasc_star:
                if s['agasc_id'] is not None:
                    star_warnings.append(
                        "Different AGASC star %s found at "
                        "given yang, zang of starcheck star %s"
                        % (agasc_star['AGASC_ID'], s['agasc_id']))
                star.mag_exp = agasc_star['MAG_ACA']
                star.mag = agasc_star['MAG_ACA']
                star.agasc_id = agasc_star['AGASC_ID']
                star.id = agasc_star['AGASC_ID']
                star.color = agasc_star['COLOR1']

            else:
                star_warnings.append("No AGASC star found for this yang, zang")

        for warning in starcheck_warnings:
            if warning.idx == star.idx:
                star_warnings.append(warning.warning)

        warnings[star.idx] = star_warnings

    return (stars, warnings)


def get_acq_data(mp_obs, stars):
    """
    For a given observation, retrieve the acquisition telemetry
    and store in the previously constructed star slot recarrays.

    :param obs: recarray as retrieved from get_needed_obsids()
    :param obs_db: recarry as retrieved from get_obs_db()
    :param stars: dict of dicts from get_acq_stars()
    """

    # retrieve the transition to NPNT from the cmd_states database if the
    # observation is in the cmd_states era
    # min_cmd_time = sqlaca.fetchone(
    # "select min(time) as time from cmds")['time']
    min_cmd_time = DateTime('2002:007:13:38:57.743').secs
    end_last_manvr_time = None
    if stars[0].kalman_tstart < mp_obs['tstart']:
        logger.warn(
            "Error, for obsid {}, mp_load_info/obidet tstart after kalman_start".format(
                mp_obs['obsid']))
        end_last_manvr_time = sqlaca.fetchone(
            """select max(tstop) as tstop from aiprops
               where tstart < %f
               and pcad_mode = 'NMAN'"""
            % (stars[0].kalman_tstart))['tstop']
    # the kalman_start and such are in all the lines of the stars recarray
    if (stars[0].kalman_tstart > min_cmd_time) and not end_last_manvr_time:
        end_last_manvr_time = sqlaca.fetchone(
            """select min(tstart) as tstart from aiprops
               where tstart < %f
               and tstart > %f
               and pcad_mode = 'NPNT'"""
            % (stars[0].kalman_tstart, mp_obs['tstart']))['tstart']
    else:
        # if before cmd_states, dig around and find the maneuver summary
        # to find the end of the maneuver before the obsid begins
        from glob import glob
        mm_files = glob(os.path.join(mp_path, mp_obs.mp_path, 'mm*.sum'))
        mm_files.extend(glob(os.path.join(mp_path + mp_obs.mp_path,
                                          'mps', 'mm*.sum')))
        mm_files.extend(glob(os.path.join(mp_path + mp_obs.mp_path,
                                          'm???:????', 'mm*.sum')))
        if not len(mm_files):
            raise ObsidError("No Maneuver Summary Found")
        import Ska.ParseCM
        mm = Ska.ParseCM.read_mm(mm_files[0])
        for m_entry in mm:
            if m_entry['tstop'] > stars[0].kalman_tstart:
                break
            end_last_manvr_time = m_entry['tstop']

    # the acq time should be sometime between the end of the last
    # maneuver and the beginning of the kalman interval
    pcad_tstart = end_last_manvr_time
    if pcad_tstart is None:
        raise ObsidError("could not determine obsid start time")
    pcad_tstop = stars[0].kalman_tstart

    acq_fields = [field + '%s' % slot for field in
                  ['AOACQID', 'AOACFCT',
                   'AOACMAG', 'AOACYAN', 'AOACZAN']
                  for slot in range(0, 8)]
    fields = ['AOACASEQ', 'COBSRQID']
    fields.extend(acq_fields)
    telem = fetch.MSIDset(fields, pcad_tstart, pcad_tstop, filter_bad=True)
    aoacaseq = telem['AOACASEQ']

    # the acquisition should also be labeled for the intended obsid,
    # so check for that
    cobsrqid = telem['COBSRQID']
    cobsid_match = np.flatnonzero(cobsrqid.vals == stars[0].obsid)
    if not len(cobsid_match):
        # If within last 4 days, throw a TooNewError that will be
        # handled by Skipping instead of Error-ing in the log
        if ((DateTime().secs - pcad_tstart) < (4 * 86400)):
            raise TooNewError(
                "obsid {} not found in cobsrqid telem in range".format(
                    stars[0].obsid))
        else:
            raise ObsidError(
                "obsid %d not found in cobsrqid telem in range"
                % stars[0].obsid)
    min_pcad_obsid = cobsrqid.times[min(cobsid_match)]
    max_pcad_obsid = cobsrqid.times[max(cobsid_match)]
    obsid_match = ((aoacaseq.times >= min_pcad_obsid)
                   & (aoacaseq.times <= max_pcad_obsid))

    gui_match = aoacaseq.vals == 'GUID'
    acq_match = aoacaseq.vals == 'AQXN'

    # find the guide transition time
    # indexes where going *to* GUID
    change_to_guid = (np.where(gui_match, 1, 0)[1:]
                      - np.where(gui_match, 1, 0)[0:-1] == 1)
    # add a false at the beginning to get the length right
    gui_trans = np.hstack([[False], change_to_guid])
    gui_idx = np.flatnonzero(gui_trans & obsid_match)
    if not len(gui_idx):
        raise ObsidError("Cannot determine guide transition time")
    if len(gui_idx) > 1:
        logger.debug("More than one guide transition")
    for poss_gui_trans in gui_idx:
        guide_time = aoacaseq.times[poss_gui_trans]
        # and the first time after that
        guide_time_plus = aoacaseq.times[poss_gui_trans + 1]
        acq_time_match = aoacaseq.times < guide_time_plus
        # and the last acquisition time before that
        acq_idx = np.flatnonzero(acq_match & obsid_match & acq_time_match)
        if len(acq_idx):
            break

    if not len(acq_idx):
        raise ObsidError("Cannot determine last ACQ time")
    acq_time = aoacaseq.times[max(acq_idx)]

    logger.debug("ACQ time %s" % DateTime(acq_time).date)

    for star in stars:
        if ((star['type'] == 'BOT') or (star['type'] == 'ACQ')):

            slot = star['slot']
            # id/noid is by catalog position...
            obc_id_msid = telem['AOACQID%s' % star.cat_pos].vals[
                (telem['AOACQID%s' % star.cat_pos].times >= guide_time)
                & (telem['AOACQID%s' % star.cat_pos].times < guide_time_plus)]

            # There was logic here for adding the made-up NOTRAK obc status
            # but that seems to be unnecessary... a star is never obc_id 'ID'
            # if it is not tracked (bad data was used in initial determination
            # that those states existed).
            star.obc_id = obc_id_msid[0].strip()

            # the fields that we want *before* the first GUID
            aoacmag = telem['AOACMAG%s' % slot].vals[
                (telem['AOACMAG%s' % slot].times >= acq_time)
                & (telem['AOACMAG%s' % slot].times < guide_time)]

            aoacyan = telem['AOACYAN%s' % slot].vals[
                (telem['AOACYAN%s' % slot].times >= acq_time)
                & (telem['AOACYAN%s' % slot].times < guide_time)]

            aoaczan = telem['AOACZAN%s' % slot].vals[
                (telem['AOACZAN%s' % slot].times >= acq_time)
                & (telem['AOACZAN%s' % slot].times < guide_time)]

            star.mag_obs = aoacmag[0]
            star.yang_obs = aoacyan[0]
            star.zang_obs = aoaczan[0]


def get_gui_data(stars, email=None):
    """
    For a given observation, retrieve the acquisition telemetry
    and store in the previously constructed star slot recarrays.

    :param obs: recarray as retrieved from get_needed_obsids()
    :param obs_db: recarry as retrieved from get_obs_db()
    :param stars: dict of dicts from get_acq_stars()
    """
    pcad_tstart = stars[0].kalman_tstart
    pcad_tstop = stars[0].kalman_tstop

    gui_fields = [field + '%s' % slot for field in [
            'AOACICC', 'AOACIDP', 'AOACIIR', 'AOACIMS',
            'AOACIQB', 'AOACISP', 'AOACQID', 'AOACFCT',
            'AOACMAG', 'AOACYAN', 'AOACZAN']
                  for slot in range(0, 8)]
    fields = ['AOACASEQ', 'COBSRQID']
    fields.extend(gui_fields)
    telem = fetch.MSIDset(fields, pcad_tstart, pcad_tstop, filter_bad=True)

    aoacaseq = telem['AOACASEQ']

    # the acquisition should also be labeled for the intended obsid,
    # so check for that
    cobsrqid = telem['COBSRQID']
    cobsid_match = np.flatnonzero(cobsrqid.vals == stars[0].obsid)
    if not len(cobsid_match):
        # If within last 4 days, throw a TooNewError that will be
        # handled by Skipping instead of Error-ing in the log
        if ((DateTime().secs - pcad_tstart) < (4 * 86400)):
            raise TooNewError(
                "obsid {} not found in cobsrqid telem in range".format(
                    stars[0].obsid))
        else:
            raise ObsidError(
                "obsid %d not found in cobsrqid telem in range"
                % stars[0].obsid)
    min_pcad_obsid = cobsrqid.times[min(cobsid_match)]
    max_pcad_obsid = cobsrqid.times[max(cobsid_match)]
    obsid_match = ((aoacaseq.times >= min_pcad_obsid)
                   & (aoacaseq.times <= max_pcad_obsid))
    for star in stars:
        if (star['type'] != 'ACQ'):
            slot = star['slot']

            refmsid = 'AOACASEQ'
            star['n_samples'] = len(telem[refmsid].times[obsid_match])
            star['not_tracking_samples'] = len(np.flatnonzero(
                telem['AOACFCT' + str(slot)].vals[obsid_match] != 'TRAK'))

            for par in ['AOACMAG', 'AOACYAN', 'AOACZAN']:
                stat_telem = telem[par + str(slot)].vals[obsid_match]
                star[par.lower() + '_min'] = np.min(stat_telem)
                star[par.lower() + '_max'] = np.max(stat_telem)
                star[par.lower() + '_mean'] = np.mean(stat_telem)
                star[par.lower() + '_median'] = np.median(stat_telem)
                star[par.lower() + '_rms'] = np.std(stat_telem)
                star[par.lower() + '_5th'] = scoreatpercentile(stat_telem, 5)
                star[par.lower() + '_95th'] = scoreatpercentile(stat_telem, 95)

            stat_map = dict(common_col='AOACICC',
                            def_pix='AOACIDP',
                            mult_star='AOACIMS',
                            ion_rad='AOACIIR',
                            quad_bound='AOACIQB',
                            sat_pix='AOACISP')

            bad_stat = np.zeros(len(telem[refmsid].times[obsid_match]),
                                dtype=bool)
            obc_bad_stat = np.zeros(len(telem[refmsid].times[obsid_match]),
                                    dtype=bool)
            obc_stat_fields = ['common_col',
                               'mult_star', 'ion_rad']
            # if observed before removing the defective pixel from the onboard filtering
            if pcad_tstart < DateTime('2013:297:11:25:52.000').secs:
                obc_stat_fields.append('def_pix')
            for imstat in stat_map.keys():
                slot_imstat = telem[stat_map[imstat] + str(slot)].vals
                obsid_slot_imstat = slot_imstat[obsid_match]
                star[imstat + '_samples'] = len(np.flatnonzero(
                        obsid_slot_imstat != 'OK '))
                bad_stat = bad_stat | (obsid_slot_imstat != 'OK ')
                if imstat in obc_stat_fields:
                    obc_bad_stat = obc_bad_stat | (obsid_slot_imstat != 'OK ')

            star['bad_status_samples'] = len(np.flatnonzero(bad_stat))
            star['obc_bad_status_samples'] = len(np.flatnonzero(obc_bad_stat))
            star['sample_interval_secs'] = np.median(
                telem[refmsid].times[obsid_match][1:]
                - telem[refmsid].times[obsid_match][0:-1])

            # warn on fids that are dropped or not tracked (more that 5%?)
            if star['type'] == 'FID':
                nt_frac = (star['not_tracking_samples'] * 1.0
                           / star['n_samples'])
                if (nt_frac > 0.05):
                    warn = (
                        "Fid in SLOT %d of " % slot
                        + "OBSID %d OBI %d " % (star.obsid, star.obi)
                        + " not tracking fraction = %.2f" % nt_frac)
                    logger.error(warn)
                    anom_email(star.obsid, star.obi, to_addr=email, mesg=warn,
                               subject='FID Trak < 95%%: Obsid %d Obi %d' % (
                            star.obsid, star.obi))


def print_debug_table(stars, warnings):

    fmt = {
        'idx': '%2d',
        'slot': '%2d',
        'type': '%s',
        'id': '%10d',
        'obc_id': '%s',
        'mag': '% 4.2f',
        'aoacmag_mean': '% 4.2f',
        'aoacyan_mean': '% 4.2f',
        'aoaczan_mean': '% 4.2f',
        'color': '% 4.2f',
        'yang': '% 8.2f',
        'zang': '% 8.2f',
        'mag_obs': '% 4.2f',
        'yang_obs': '% 8.2f',
        'zang_obs': '% 8.2f',
        'ap_date': '%s',
        'revision': '%s',
        'd_mag': '% 4.2f',
        'd_yang': '% 8.2f',
        'd_zang': '% 8.2f',
        }

    fields = ['idx', 'slot', 'type', 'id', 'obc_id',
              'mag', 'mag_obs', 'aoacmag_mean',
              'yang', 'yang_obs', 'aoacyan_mean',
              'zang', 'zang_obs', 'aoaczan_mean']
    logger.debug("\t".join(fields))
    for star in stars:
        slot_line = []
        for field in fields:
            try:
                slot_line.append(fmt[field] % getattr(star, field))
                slot_line.append("\t")
            except TypeError:
                slot_line.append("N/A\t")
        logger.debug("".join(slot_line))

    for idx in stars['idx']:
        for warn in warnings[idx]:
                logger.debug("idx: % 2d %s" % (idx, warn))


def update_db(stars, warnings, dbh):
    """
    Delete guide stats entries for the given obsid and insert the new entries.
    """

    acqs = stars[(stars['type'] == 'BOT')
                 | (stars['type'] == 'ACQ')][acq_cols]
    if len(acqs):
        dbh.execute("delete from %s where obsid = %d and obi = %d"
                    % (data_table['acq'], acqs[0]['obsid'], acqs[0]['obi']))
        dbh.execute("delete from %s where obsid = %d and obi = %d"
                    % (warning_table['acq'], acqs[0]['obsid'], acqs[0]['obi']))
        logger.info("Updating ACQ stars for obsid = %d" % acqs[0]['obsid'])
        logger.debug("Inserting acq stars")
        for acq in acqs:
            if len(warnings[acq['idx']]):
                for warn in warnings[acq['idx']]:
                    dbh.insert(dict(obsid=acq['obsid'],
                                    obi=acq['obi'],
                                    slot=acq['slot'],
                                    details=warn), warning_table['acq'])
            dbh.insert(acq, data_table['acq'])
        logger.debug("acq inserts complete")
    trak = stars[(stars['type'] == 'FID')
                 | (stars['type'] == 'GUI')
                 | (stars['type'] == 'BOT')][gui_cols]
    if len(trak):
        dbh.execute("delete from %s where obsid = %d and obi = %d"
                    % (data_table['gui'], trak[0]['obsid'], trak[0]['obi']))
        dbh.execute("delete from %s where obsid = %d and obi = %d"
                    % (warning_table['gui'], trak[0]['obsid'], trak[0]['obi']))
        logger.info("Updating GUI stars for obsid = %d" % trak[0]['obsid'])
        logger.debug("Inserting gui stars")
        for star in trak:
            if len(warnings[star['idx']]):
                for warn in warnings[star['idx']]:
                    dbh.insert(dict(obsid=star['obsid'],
                                    obi=star['obi'],
                                    slot=star['slot'],
                                    warning=warn), warning_table['gui'])
            dbh.insert(star, data_table['gui'])
        logger.debug("gui inserts complete")

def anom_email(obsid, obi, to_addr, mesg, subject):
    msg = MIMEText(mesg)
    msg['From'] = 'aca@head.cfa.harvard.edu'
    msg['Subject'] = subject
    msg['To'] = to_addr
    s = smtplib.SMTP('head.cfa.harvard.edu')
    s.sendmail('aca@head.cfa.harvard.edu', [to_addr], msg.as_string())
    s.quit()


def check_acq_id_count(stars, limit=3, email=None):
    """
    For each obsid, count the number of acquired stars and
    throw a warning if at or under limit.

    :param stars: stars structure from get_acq_stars()
    :param limit: threshold for warnings
    :param email: email address for warnings
    """
    id_count = np.count_nonzero(stars['obc_id'] == 'ID')
    if id_count <= limit:
        anom_text = "Warning: Obsid {} shows only {} identified ACQ stars".format(
            stars[0]['obsid'], id_count)
        if email:
            anom_email(stars[0]['obsid'], stars[0]['obi'], to_addr=email, mesg=anom_text,
                       subject='Only {} ID\'d acq stars for Obsid {}'.format(
                    id_count,
                    stars[0]['obsid']))
        logger.error(anom_text)


def get_acq_deltas(stars, email=None):
    """
    For the obsid, and the stars structure given, calculate
    the differences between the observed and expected mag, yang, and zang.

    :param obs: obs recarray from get_needed_obsids()
    :param stars: stars structure created in get_acq_stars() and
    populated with telem from get_pcad_data()
    """

    for star in stars:
        # if star has magnitude data, consider it 'tracked'
        if ((star.mag_obs < 13.9) &
            ((star['type'] == 'BOT') | (star['type'] == 'ACQ'))):
            # calculate mean offset of all of the other tracked stars
            y_sum = 0
            z_sum = 0
            number = 0
            ok_mask = (((stars['type'] == 'ACQ') | (stars['type'] == 'BOT'))
                       & (stars['mag_obs'] < 13.9) & (stars['obc_id'] == 'ID'))
            for o_star in stars[ok_mask]:
                if o_star.slot != star.slot:
                    y_sum += o_star.yang_obs - o_star.yang
                    z_sum += o_star.zang_obs - o_star.zang
                    number += 1

            # using the mean offset, calculated the differences
            # in expected/observed positions
            if number:
                y_off = y_sum / number
                z_off = z_sum / number
                star.d_yang = star.yang_obs - (star.yang + y_off)
                star.y_offset = y_off
                star.d_zang = star.zang_obs - (star.zang + z_off)
                star.z_offset = z_off
            else:
                raise ObsidError("No stars to calculate mean offset")
            star.d_mag = star.mag_obs - star.mag
            d_rad = (star.d_yang ** 2 + star.d_zang ** 2) ** .5
            if d_rad >= acq_anom_radius:
                anom_text = ("Large Deviation from Expected ACQ Star Position "
                             + "in %d (obi %d)\n" % (star.obsid, star.obi)
                             + "\tSlot %d Expected (Y-Pos,Z-Pos) = (%.1f, %.1f) \n"
                             % (star.slot, star.yang, star.zang)
                             + "\tSlot %d Observed (Y-Pos,Z-Pos) = (%.1f, %.1f) \n"
                             % (star.slot, star.yang_obs, star.zang_obs))

                known_anoms = Ska.Table.read_ascii_table(anom_list)
                logger.debug(
                    "Reading list of known acq anoms from %s" % anom_list)
                anoms = set((x['obsid'], x['obi'], x['slot'])
                            for x in known_anoms)
                curr_anom = (star['obsid'], star['obi'], star['slot'])
                # if we've already seen it, don't "warn" and send email
                if curr_anom in anoms:
                    logger.info(anom_text)
                else:
                    if email:
                        anom_email(star.obsid, star.obi, to_addr=email, mesg=anom_text,
                                       subject='Acq Anomaly: Obsid %d Obi %d' % (
                                star.obsid, star.obi))
                        logger.error(anom_text)

def update_obi(obs, dbh, dryrun=False, email=None):
    obs_db = get_obs_db(obs.obsid, obs.obi, obs.tstart)
    stars, warnings = get_stars(obs_db, obs, dbh.dbi)
    get_acq_data(obs, stars)
    get_acq_deltas(stars, email)
    check_acq_id_count(stars, email=email)
    get_gui_data(stars, email)
    print_debug_table(stars, warnings)
    if not dryrun:
        update_db(stars, warnings, dbh=dbh)


def main():

    (opt,args) = get_options()

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    if opt.verbose == 2:
        ch.setLevel(logging.DEBUG)
    if opt.verbose == 0:
        ch.setLevel(logging.WARN)
    if not len(logger.handlers):
        logger.addHandler(ch)

    nowdate=time.ctime()
    logger.info("---------- star stats DB update at %s ----------" % (nowdate))
    # make a sqlite db if in that mode and one doesn't exist
    if (opt.dbi == 'sqlite' and
        (not os.path.exists(opt.server) or os.stat(opt.server).st_size == 0)):
        logger.info("Creating db from create_tables.sqlite")
        db_init_cmds = file("create_tables.sqlite").read()
        db = Ska.DBI.DBI(dbi='sqlite', server=opt.server)
        db.execute(db_init_cmds)
        del db
    dbh = Ska.DBI.DBI(dbi=opt.dbi, server=opt.server, user=opt.user, database=opt.database)
    logger.debug("connecting to db (dbi=%s, server=%s, user=%s, database=%s)" %
                (dbh.dbi, dbh.server, dbh.user, dbh.database))
    okmissing = Ska.Table.read_ascii_table(opt.missing_list)
    logger.debug("Reading list of expected missing from %s" % opt.missing_list)
    okskip = set((x['obsid'], x['obi']) for x in okmissing)
    obsdata = get_needed_obsids(requested_obsid=opt.obsid,
                                missing_set=okskip)

    for obs in obsdata:
        logger.debug("Processing %d %d" % (obs.obsid, obs.obi))
        try:
            update_obi(obs, dbh=dbh, dryrun=opt.dryrun, email=opt.email)
        except TooNewError as detail:
            logger.debug("Skipping Too New obsid %d obi %d: %s" % (obs.obsid,
                                                                  obs.obi,
                                                                  detail))
        except (TooOldError, WeirdObsidError) as detail:
            okskip = okskip | set([(obs.obsid, obs.obi)])
            logger.debug("Skipping obsid %d obi %d: %s" % (obs.obsid,
                                                           obs.obi,
                                                           detail))

        except (MissingDataError, ObsidError, ExceptionPexpect) as detail:
            logger.warn("Error processing %d:%d obidet date: '%s' : %s"
                        % (obs.obsid,
                           obs.obi,
                           obs.last_ap_date,
                           detail))
    if opt.update_missing:
        f = open(opt.missing_list, 'w')
        f.write("obsid,obi,date\n")
        for obsid, obi in okskip:
            f.write("%s,%s,'%s'\n" % (obsid, obi))
        f.close()

    logger.removeHandler(ch)

if __name__ == '__main__':
    main()
