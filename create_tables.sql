create table acq_stats_data
(obsid		int		not null,
obi		int		not null,
tstart		float,
tstop		float,
slot		int		not null,
idx int,
cat_pos int,
type	    varchar(3) not null,
agasc_id	int		not null,
obc_id      varchar(6) not null,
yang		float,
zang		float,
mag		    float,
color float,
halfw float,
mag_obs		float,
yang_obs	float,
zang_obs	float,
d_mag       float,
d_yang      float,
d_zang      float,
y_offset float,
z_offset float,
ap_date      datetime,
revision varchar(10))
go
sp_primarykey acq_stats_data, obsid, obi, slot
go

create table acq_stats_warnings
(obsid           int     not null,
obi		int   not null,
slot            int  not null,
details		varchar(100))
go





create table trak_stats_data
(obsid		int		not null,
obi		int		not null,
slot		int		not null,
idx	int,
cat_pos	int,
id int,
type varchar(4),
color float,
cyan_exp float,
czan_exp float,
mag_exp float,
kalman_datestart varchar(21),
kalman_datestop varchar(21),
kalman_tstart   float,
kalman_tstop    float,
aoacyan_min	float,
aoacyan_mean	float,
aoacyan_max	float,
aoacyan_rms	float,
aoacyan_median	float,
aoacyan_5th	float,
aoacyan_95th	float,
aoaczan_min	float,
aoaczan_mean	float,
aoaczan_max	float,
aoaczan_rms	float,
aoaczan_median	float,
aoaczan_5th	float,
aoaczan_95th	float,
aoacmag_min	float,
aoacmag_mean	float,
aoacmag_max	float,
aoacmag_rms	float,
aoacmag_median	float,
aoacmag_5th	float,
aoacmag_95th	float,
n_samples	int,
not_tracking_samples	int,
bad_status_samples	int,
obc_bad_status_samples  int,
common_col_samples	int,
sat_pix_samples		int,
def_pix_samples   	int,
quad_bound_samples	int,
ion_rad_samples   	int,
mult_star_samples	int,
sample_interval_secs	float,
ap_date      datetime,
revision varchar(10))
go
sp_primarykey trak_stats_data, obsid, obi, slot
go

create table trak_slot_warnings
(obsid		int		not null,
obi		int		not null,
slot		int		not null,
idx	int,
id int,
type varchar(4),
warning varchar(100),
)
go


