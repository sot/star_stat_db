drop view guide_stats_view;
create view guide_stats_view as
select
gsd.obsid as obsid,
gsd.obi as obi,
gsd.slot as slot,
gsd.kalman_tstart as kalman_tstart,
gsd.kalman_tstop as kalman_tstop,
sc.yang as yang_exp,
sc.zang as zang_exp,
sc.mag as mag_exp,
sc.type as type,
sc.id as id,
aoacyan_min as yang_obs_min,
aoacyan_mean as yang_obs_mean,
aoacyan_max as yang_obs_max,
aoacyan_rms as yang_obs_rms,
aoacyan_median as yang_obs_50th,
aoacyan_5th as yang_obs_5th,
aoacyan_95th as yang_obs_95th,
aoaczan_min as zang_obs_min,
aoaczan_mean as zang_obs_mean,
aoaczan_max as zang_obs_max,
aoaczan_rms as zang_obs_rms,
aoaczan_median as zang_obs_50th,
aoaczan_5th as zang_obs_5th,
aoaczan_95th as zang_obs_95th,
aoacmag_min as mag_obs_min,
aoacmag_mean as mag_obs_mean,
aoacmag_max as mag_obs_max,
aoacmag_rms as mag_obs_rms,
aoacmag_median as mag_obs_50th,
aoacmag_5th as mag_obs_5th,
aoacmag_95th as mag_obs_95th,
color,
n_samples,
not_tracking_samples,
obc_bad_status_samples,
bad_status_samples,
common_col_samples,
sat_pix_samples,
def_pix_samples,
quad_bound_samples,
ion_rad_samples,
mult_stars_samples,
sample_interval_secs,
(not_tracking_samples*1.0/n_samples)*100  as percent_not_tracking,
(obc_bad_status_samples*1.0/n_samples)*100 as percent_obc_bad_status,
(bad_status_samples*1.0/n_samples)*100  as percent_bad_status,
(common_col_samples*1.0/n_samples)*100  as percent_common_col,
(sat_pix_samples*1.0/n_samples)*100  as percent_sat_pix,
(def_pix_samples*1.0/n_samples)*100  as percent_def_pix,
(quad_bound_samples*1.0/n_samples)*100  as percent_quad_bound,
(ion_rad_samples*1.0/n_samples)*100  as percent_ion_rad,
(mult_stars_samples*1.0/n_samples)*100  as percent_mult_stars,
gsd.ap_date as ap_date,
gsd.tool_cvs_rev as tool_cvs_rev
from
guide_stats_data as gsd,
starcheck_catalog as sc 
where 
gsd.obsid = sc.obsid
and gsd.obi = sc.obi
and gsd.slot = sc.slot
and ( sc.type = 'BOT' or sc.type = 'GUI' or sc.type = 'FID' );

drop view guide_status_percent;
create view guide_status_percent as
select
gsd.obsid as obsid,
gsd.obi as obi,
gsd.slot as slot,
sc.idx as idx,
sc.type as type,
sc.id as id,
str(((not_tracking_samples*1.0/n_samples)*100),6,2) as no_trak,
str(((obc_bad_status_samples*1.0/n_samples)*100),6,2) as obc_bad_stat,
str(((bad_status_samples*1.0/n_samples)*100),6,2) as bad_stat,
str(((common_col_samples*1.0/n_samples)*100),6,2) as cc,
str(((sat_pix_samples*1.0/n_samples)*100),6,2) as sp,
str(((def_pix_samples*1.0/n_samples)*100),6,2) as dp,
str(((quad_bound_samples*1.0/n_samples)*100),6,2) as qb,
str(((ion_rad_samples*1.0/n_samples)*100),6,2) as ir,
str(((mult_stars_samples*1.0/n_samples)*100),6,2) as ms
from
guide_stats_data as gsd,
starcheck_catalog as sc 
where 
gsd.obsid = sc.obsid
and gsd.obi = sc.obi
and gsd.slot = sc.slot
and ( sc.type = 'BOT' or sc.type = 'GUI' or sc.type = 'FID' );

drop view guide_status_counts;
create view guide_status_counts as
select
gsd.obsid as obsid,
gsd.obi as obi,
gsd.slot as slot,
sc.type as type,
sc.id as id,
n_samples,
not_tracking_samples,
obc_bad_status_samples,
bad_status_samples,
common_col_samples,
sat_pix_samples,
def_pix_samples,
quad_bound_samples,
ion_rad_samples,
mult_stars_samples
from
guide_stats_data as gsd,
starcheck_catalog as sc 
where 
gsd.obsid = sc.obsid
and gsd.obi = sc.obi
and gsd.slot = sc.slot
and ( sc.type = 'BOT' or sc.type = 'GUI' or sc.type = 'FID' );

drop view guide_stats_ranges;
create view guide_stats_ranges as
select
gsd.obsid as obsid,
gsd.obi as obi,
gsd.slot as slot,
sc.yang as yang_exp,
sc.zang as zang_exp,
sc.mag as mag_exp,
sc.type as type,
sc.id as id,
aoacyan_min as yang_obs_min,
aoacyan_mean as yang_obs_mean,
aoacyan_max as yang_obs_max,
aoacyan_rms as yang_obs_rms,
aoacyan_median as yang_obs_50th,
aoacyan_5th as yang_obs_5th,
aoacyan_95th as yang_obs_95th,
aoaczan_min as zang_obs_min,
aoaczan_mean as zang_obs_mean,
aoaczan_max as zang_obs_max,
aoaczan_rms as zang_obs_rms,
aoaczan_median as zang_obs_50th,
aoaczan_5th as zang_obs_5th,
aoaczan_95th as zang_obs_95th,
aoacmag_min as mag_obs_min,
aoacmag_mean as mag_obs_mean,
aoacmag_max as mag_obs_max,
aoacmag_rms as mag_obs_rms,
aoacmag_median as mag_obs_50th,
aoacmag_5th as mag_obs_5th,
aoacmag_95th as mag_obs_95th,
color
from
guide_stats_data as gsd,
starcheck_catalog as sc 
where 
gsd.obsid = sc.obsid
and gsd.obi = sc.obi
and gsd.slot = sc.slot
and ( sc.type = 'BOT' or sc.type = 'GUI' or sc.type = 'FID' );


