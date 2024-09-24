#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* prototypes */
int readmat(char *fname, float ***A, int *rows, int *cols);
int gentrap(float A, float y_max, float dydx_max, float dx, float *y, float *x_ramp, float *x_plat);

int gengrads(int id, float dt, float gmax, float smax, /* inputs: trajectory id #, raster time (s), max gradient amp (G/cm), max slew rate (G/cm/s) */
		float ***gx, float ***gy, float ***gz, /* output gradient waveforms: gx, gy, gz (G/cm, all of size [nshots x n_all]) */
		int *n_shots, int *n_all, int *n_wav, int *n_off /* output sizes: # of shots, number of points in gradient waveforms, number of points in trajectory, acq offset */
		)
{
/* gengrads() by David Frey (2024):
 * 
 * Generates readout gradients for arbitrary kspace trajectory from file
 *
 * kspace trajectory should be contained in 3 files: sparky_trajectories/ {sparky#####.kx, sparky#####.ky and sparky#####.kz}
 * each file should be a matrix of size [n_wav x n_shots] for corresponding kx, ky, or kz waveforms
 * NOTE: be sure to remove extra space at end of each row; each value should be tab separated
 *
 */

	/* define constants */
	float gam = 4258; /* gyromagnetic ratio (Hz/G) */

	/* declare waveforms */
	float **kx, **ky, **kz; /* kspace trajectory (2D array structured like a matrix - not array of arrays) */
	float **gx_wav, **gy_wav, **gz_wav; /* gradients through kspace */

	/* declare variables */
	float *gx0, *gy0, *gz0, *gxend, *gyend, *gzend; /* start & end wav gradient amplitudes */
	float g0_max, gend_max;	/* max start & end wav gradient amplitudes */
	float *dkx0, *dky0, *dkz0, *dkxend, *dkyend, *dkzend; /* start & end kspace after ramps */
	float dk0_max, dkend_max; /* max start & end kspace after ramps */
	float t_rup, t_rdn; /* gradient ramp times */
	float a_pwd_max, t_pwd_ramp, t_pwd_plat; /* prewinder trapezoid times */
	float a_rwd_max, t_rwd_ramp, t_rwd_plat; /* rewinder trapezoid times */
	float ax_rwd, ay_rwd, az_rwd, ax_pwd, ay_pwd, az_pwd;
	int n_pwd_ramp, n_pwd_plat, n_rup, n_rdn, n_rwd_ramp, n_rwd_plat; /* number of samples in ramps/traps */
	int n0, n, shot; /* loop indicies */
	int tmp_n_wav, tmp_n_shots; /* matrix sizes */
	char fname[100]; /* file name buffer */

	/* write out the id number to file */
	FILE *fID_sparky_trajid = fopen("sparky_trajid.txt","w");
	fprintf(fID_sparkytrajid, "%05d", id);
	fclose(fID_sparky_trajid);

	/* by default - set gradients to 0 to acquire FID only */
	if (!id) {
		fprintf(stderr, "gengrads(): no id passed, setting gradients to 0\n");

		/* FID of length 4ms, with some extra initial/post buffer */
		*n_shots = 10;
		*n_wav = 1000;
		*n_off = 100;
		*n_all = 1200;
		
		/* allocate memory for the full gradient waveforms */
		*gx = (float **)malloc(*n_shots * sizeof(float *));
		*gy = (float **)malloc(*n_shots * sizeof(float *));
		*gz = (float **)malloc(*n_shots * sizeof(float *));
		for (shot = 0; shot < *n_shots; shot++) {
			(*gx)[shot] = (float *)malloc(*n_all * sizeof(float));
			(*gy)[shot] = (float *)malloc(*n_all * sizeof(float));
			(*gz)[shot] = (float *)malloc(*n_all * sizeof(float));
		}

		/* set gradients to 0 */
		for (shot = 0; shot < *n_shots; shot++) {
			for (n = 0; n < *n_all; n++) {
				(*gx)[shot][n] = 0;
				(*gy)[shot][n] = 0;
				(*gz)[shot][n] = 0;
			}
		}

		return 1;
	}

	/* read in the kx waveform (n_wav x n_shots) */
	sprintf(fname, "./sparky_trajectories/sparky%05d.kx", id);
	if (!readmat(fname, &kx, &tmp_n_wav, &tmp_n_shots)) {
		fprintf(stderr, "gengrads(): failed to read matrix from file %s\n", fname);
		return 0;
	}
	*n_wav = tmp_n_wav;
	*n_shots = tmp_n_shots;

	/* read in the ky waveform (n_wav x n_shots) */
	sprintf(fname, "./sparky_trajectories/sparky%05d.ky", id);
	if (!readmat(fname, &ky, &tmp_n_wav, &tmp_n_shots)) {
		fprintf(stderr, "gengrads(): failed to read matrix from file %s\n", fname);
		return 0;
	}
	if (tmp_n_wav != *n_wav || tmp_n_shots != *n_shots) {
		fprintf(stderr, "gengrads(): matrix size discrepency in %s: (%d x %d) vs kx file: (%d x %d)\n", 
				fname, tmp_n_wav, tmp_n_shots, *n_wav, *n_shots);
		return 0;
	}

	/* read in the kz waveform (n_wav x n_shots) */
	sprintf(fname, "./sparky_trajectories/sparky%05d.kz", id);
	if (!readmat(fname, &kz, &tmp_n_wav, &tmp_n_shots)) {
		fprintf(stderr, "gengrads(): failed to read matrix from file %s\n", fname);
		return 0;
	}
	if (tmp_n_wav != *n_wav || tmp_n_shots != *n_shots) {
		fprintf(stderr, "gengrads(): matrix size discrepency in %s: (%d x %d) vs kx file: (%d x %d)\n", 
				fname, tmp_n_wav, tmp_n_shots, *n_wav, *n_shots);
		return 0;
	}

	/* allocate memory for the gradient waveforms without ramps */
	gx_wav = (float **)malloc(*n_shots * sizeof(float *));
	gy_wav = (float **)malloc(*n_shots * sizeof(float *));
	gz_wav = (float **)malloc(*n_shots * sizeof(float *));
	for (shot = 0; shot < *n_shots; shot++) {
		gx_wav[shot] = (float *)malloc((*n_wav-1) * sizeof(float));
		gy_wav[shot] = (float *)malloc((*n_wav-1) * sizeof(float));
		gz_wav[shot] = (float *)malloc((*n_wav-1) * sizeof(float));
	}
	
	/* allocate memory for the initial/final values to store */
	gx0 = (float *)malloc(*n_shots * sizeof(float));
	gy0 = (float *)malloc(*n_shots * sizeof(float));
	gz0 = (float *)malloc(*n_shots * sizeof(float));
	gxend = (float *)malloc(*n_shots * sizeof(float));
	gyend = (float *)malloc(*n_shots * sizeof(float));
	gzend = (float *)malloc(*n_shots * sizeof(float));
	dkx0 = (float *)malloc(*n_shots * sizeof(float));
	dky0 = (float *)malloc(*n_shots * sizeof(float));
	dkz0 = (float *)malloc(*n_shots * sizeof(float));
	dkxend = (float *)malloc(*n_shots * sizeof(float));
	dkyend = (float *)malloc(*n_shots * sizeof(float));
	dkzend = (float *)malloc(*n_shots * sizeof(float));

	/* calculate the gradients, determine max start & end gradients */
	g0_max = 0;
	gend_max = 0;
	for (shot = 0; shot < *n_shots; shot++) {
		for (n = 0; n < *n_wav-1; n++) {

			/* calculate x gradient */
			gx_wav[shot][n] = (kx[n+1][shot] - kx[n][shot]) / (gam*dt); 
			if (fabs(gx_wav[shot][n]) > gmax) {
				fprintf(stderr, "gengrads(): kx exceeds gradient amp limit (%f G/cm) at point %d, shot %d\n", gmax, n, shot);
				return 0;
			}
			if (n > 1 && fabs(gx_wav[shot][n] - gx_wav[shot][n-1])/dt > smax) {
				fprintf(stderr, "gengrads(): kx exceeds slew limit (%f G/cm/s) at point %d, shot %d\n", smax, n, shot);
				return 0;
			}

			/* calculate y gradient */
			gy_wav[shot][n] = (ky[n+1][shot] - ky[n][shot]) / (gam*dt);
			if (fabs(gy_wav[shot][n]) > gmax) {
				fprintf(stderr, "gengrads(): ky exceeds gradient amp limit (%f G/cm) at point %d, shot %d\n", gmax, n, shot);
				return 0;
			}
			if (n > 1 && fabs(gy_wav[shot][n] - gy_wav[shot][n-1])/dt > smax) {
				fprintf(stderr, "gengrads(): ky exceeds slew limit (%f G/cm/s) at point %d, shot %d\n", smax, n, shot);
				return 0;
			}

			/* calculate z gradient */
			gz_wav[shot][n] = (kz[shot][n+1] - kz[shot][n]) / (gam*dt);
			if (fabs(gz_wav[shot][n]) > gmax) {
				fprintf(stderr, "gengrads(): kz exceeds gradient amp limit (%f G/cm) at point %d, shot %d\n", gmax, n, shot);
				return 0;
			}
			if (n > 1 && fabs(gz_wav[shot][n] - gz_wav[shot][n-1])/dt > smax) {
				fprintf(stderr, "gengrads(): kz exceeds slew limit (%f G/cm/s) at point %d, shot %d\n", smax, n, shot);
				return 0;
			}

		}

		/* get initial gradients and store max */
		gx0[shot] = gx_wav[shot][0];
		g0_max = fmax(fabs(g0_max), fabs(gx0[shot]));
		gy0[shot] = gy_wav[shot][0];
		g0_max = fmax(fabs(g0_max), fabs(gy0[shot]));
		gz0[shot] = gz_wav[shot][0];
		g0_max = fmax(fabs(g0_max), fabs(gy0[shot]));
		
		/* get final gradients and store max */
		gxend[shot] = gx_wav[shot][*n_wav-2];
		gend_max = fmax(fabs(gend_max), fabs(gxend[shot]));
		gyend[shot] = gy_wav[shot][*n_wav-2];
		gend_max = fmax(fabs(gend_max), fabs(gyend[shot]));
		gzend[shot] = gz_wav[shot][*n_wav-2];
		gend_max = fmax(fabs(gend_max), fabs(gyend[shot]));

	}

	/* determine ramp times */
	t_rup = dt * ceil(g0_max / smax / dt);
	t_rdn = dt * ceil(gend_max / smax / dt);

	/* loop through shots and determine max pre/rewinder area */
	dk0_max = 0;
	dkend_max = 0;
	for (shot = 0; shot < *n_shots; shot++) {
		
		/* calculate initial kspace displacement, including accumulation from ramp-up */
		dkx0[shot] = kx[0][shot] - 0.5*gx0[shot]*gam*t_rup;
		dk0_max = fmax(fabs(dk0_max), fabs(dkx0[shot]));
		dky0[shot] = ky[0][shot] - 0.5*gy0[shot]*gam*t_rup;
		dk0_max = fmax(fabs(dk0_max), fabs(dky0[shot]));
		dkz0[shot] = kz[0][shot] - 0.5*gz0[shot]*gam*t_rup;
		dk0_max = fmax(fabs(dk0_max), fabs(dkz0[shot]));
		
		/* calculate final kspace displacement, including accumulation from ramp-down */
		dkxend[shot] = kx[*n_wav-1][shot] + 0.5*gxend[shot]*gam*t_rup;
		dkend_max = fmax(fabs(dkend_max), fabs(dkxend[shot]));
		dkyend[shot] = ky[*n_wav-1][shot] + 0.5*gyend[shot]*gam*t_rup;
		dkend_max = fmax(fabs(dkend_max), fabs(dkyend[shot]));
		dkzend[shot] = kz[*n_wav-1][shot] + 0.5*gzend[shot]*gam*t_rup;
		dkend_max = fmax(fabs(dkend_max), fabs(dkzend[shot]));

	}

	/* determine the prewinder and rewinder timings/max amplitudes */
	gentrap(dk0_max/gam, gmax, smax, dt, &a_pwd_max, &t_pwd_ramp, &t_pwd_plat);
	gentrap(dkend_max/gam, gmax, smax, dt, &a_rwd_max, &t_rwd_ramp, &t_rwd_plat);

	/* convert times to number of samples */
	*n_all = 0;
	*n_off = 0;
	n_pwd_ramp = round(t_pwd_ramp / dt);
	*n_all += 2*n_pwd_ramp;
	*n_off += 2*n_pwd_ramp;
	n_pwd_plat = round(t_pwd_plat / dt);
	*n_all += n_pwd_plat;
	*n_off += 2*n_pwd_plat;
	n_rup = round(t_rup / dt);
	*n_all += n_rup;
	*n_off += n_rup;
	*n_all += *n_wav-1; /* kspace gradients */
	n_rdn = round(t_rdn / dt);
	*n_all += n_rdn;
	n_rwd_ramp = round(t_rwd_ramp / dt);
	*n_all += 2*n_rwd_ramp;
	n_rwd_plat = round(t_rwd_plat / dt);
	*n_all += n_rwd_plat;
	
	/* allocate memory for the full gradient waveforms */
	*gx = (float **)malloc(*n_shots * sizeof(float *));
	*gy = (float **)malloc(*n_shots * sizeof(float *));
	*gz = (float **)malloc(*n_shots * sizeof(float *));
	for (shot = 0; shot < *n_shots; shot++) {
		(*gx)[shot] = (float *)malloc(*n_all * sizeof(float));
		(*gy)[shot] = (float *)malloc(*n_all * sizeof(float));
		(*gz)[shot] = (float *)malloc(*n_all * sizeof(float));
	}

	/* append ramps and traps to each waveform */
	for (shot = 0; shot < *n_shots; shot++) {

		/* calculate prewinder amplitudes for current shot */
		ax_pwd = dkx0[shot]/dk0_max * a_pwd_max;
		ay_pwd = dky0[shot]/dk0_max * a_pwd_max;
		az_pwd = dkz0[shot]/dk0_max * a_pwd_max;

		/* calculate rewinder amplitudes for current shot */
		ax_rwd = -dkxend[shot]/dkend_max * a_rwd_max;
		ay_rwd = -dkyend[shot]/dkend_max * a_rwd_max;
		az_rwd = -dkzend[shot]/dkend_max * a_rwd_max;

		/* initialize n0 */
		n0 = 0;

		/* prewinder ramp-up */
		for (n = n0; n < n0+n_pwd_ramp; n++) {
			(*gx)[shot][n] = ax_pwd * (float)(n-n0)/n_pwd_ramp;
			(*gy)[shot][n] = ay_pwd * (float)(n-n0)/n_pwd_ramp;
			(*gz)[shot][n] = az_pwd * (float)(n-n0)/n_pwd_ramp;
		}
		n0 = n;
		
		/* prewinder plateau */
		for (n = n0; n < n0+n_pwd_plat; n++) {
			(*gx)[shot][n] = ax_pwd;
			(*gy)[shot][n] = ay_pwd;
			(*gz)[shot][n] = az_pwd;
		}
		n0 = n;
		
		/* prewinder ramp-down */
		for (n = n0; n < n0+n_pwd_ramp; n++) {
			(*gx)[shot][n] = ax_pwd * (1 - (float)(n-n0)/n_pwd_ramp);
			(*gy)[shot][n] = ay_pwd * (1 - (float)(n-n0)/n_pwd_ramp);
			(*gz)[shot][n] = az_pwd * (1 - (float)(n-n0)/n_pwd_ramp);
		}
		n0 = n;

		/* kspace gradients ramp-up */
		for (n = n0; n < n0+n_rup; n++) {
			(*gx)[shot][n] = gx0[shot] * (float)(n-n0)/n_rup;
			(*gy)[shot][n] = gy0[shot] * (float)(n-n0)/n_rup;
			(*gz)[shot][n] = gz0[shot] * (float)(n-n0)/n_rup;
		}
		n0 = n;

		/* kspace gradients */
		for (n = n0; n < n0+ *n_wav-1; n++) {
			(*gx)[shot][n] = gx_wav[shot][n-n0];
			(*gy)[shot][n] = gy_wav[shot][n-n0];
			(*gz)[shot][n] = gz_wav[shot][n-n0];
		}
		n0 = n;
		
		/* kspace gradients ramp-down */
		for (n = n0; n < n0+n_rdn; n++) {
			(*gx)[shot][n] = gxend[shot] * (1 - (float)(n-n0)/n_rdn);
			(*gy)[shot][n] = gyend[shot] * (1 - (float)(n-n0)/n_rdn);
			(*gz)[shot][n] = gzend[shot] * (1 - (float)(n-n0)/n_rdn);
		}
		n0 = n;

		/* rewinder ramp-up */
		for (n = n0; n < n0+n_rwd_ramp; n++) {
			(*gx)[shot][n] = ax_rwd * (float)(n-n0)/n_rwd_ramp;
			(*gy)[shot][n] = ay_rwd * (float)(n-n0)/n_rwd_ramp;
			(*gz)[shot][n] = az_rwd * (float)(n-n0)/n_rwd_ramp;
		}
		n0 = n;

		/* rewinder plateau */
		for (n = n0; n < n0+n_rwd_plat; n++) {
			(*gx)[shot][n] = ax_rwd;
			(*gy)[shot][n] = ay_rwd;
			(*gz)[shot][n] = az_rwd;
		}
		n0 = n;
		
		/* rewinder ramp-down */
		for (n = n0; n < n0+n_rwd_ramp; n++) {
			(*gx)[shot][n] = ax_rwd * (1 - (float)(n-n0)/n_rwd_ramp);
			(*gy)[shot][n] = ay_rwd * (1 - (float)(n-n0)/n_rwd_ramp);
			(*gz)[shot][n] = az_rwd * (1 - (float)(n-n0)/n_rwd_ramp);
		}
		n0 = n;
	
	}

	return 1;
}

int readmat(char *fname, float ***A, int *rows, int *cols) {

	/* open the file or return error otherwise */
	FILE *fID = fopen(fname, "r");
	if (!fID) {
		fprintf(stderr, "readmat(): could not open file %s\n", fname);
		return 0;
	}

	char ch;
	int current_cols = 0;
	int i, j;
	*rows = 0;
	*cols = 0;

	/* first pass: determine the size of the matrix (rows and columns) */
	while ((ch = fgetc(fID)) != EOF) {
		if (ch == '\t') {
			current_cols++;
		} else if (ch == '\n') {
			(*rows)++;
			if (*cols == 0) {
				*cols = current_cols + 1;
			}
			current_cols = 0;
		}
	}
	rewind(fID);

	/* allocate memory */
	*A = (float **)malloc((*rows+1) * sizeof(float *));
	for (i = 0; i < *rows+1; i++) {
		(*A)[i] = (float *)malloc((*cols+1) * sizeof(float));
	}

	/* second pass: read the data */
	for (i = 0; i < *rows; i++) {
		for (j = 0; j < *cols; j++) {
			fscanf(fID, "%f", &(*A)[i][j]);
			if (j < *cols - 1) {
				fgetc(fID);
			}
		}
		fgetc(fID);
	}

	fclose(fID);
	return 1;
}

int gentrap(float A, float y_max, float dydx_max, float dx, float *y, float *x_ramp, float *x_plat) {

	/* determine maximum ramp time at full slew */
	float x_ramp_max = y_max / dydx_max;
	x_ramp_max = dx * ceil(x_ramp_max / dx); /* round to nearest sampling interval */

	if (A > x_ramp_max*y_max) { /* plateau needed */
		*x_ramp = x_ramp_max;
		*x_plat = (A - (x_ramp_max*y_max)) / y_max;
		*x_plat = dx * ceil(*x_plat / dx);
		*y = A / (*x_plat + *x_ramp); /* account for rounding */	
	}
	else { /* no plateau needed - only triangle */
		*x_plat = 0;
		*x_ramp = sqrt(A / dydx_max);
		*x_ramp = dx * ceil(*x_ramp / dx);
		*y = A / *x_ramp;
	}

	return 1;
}
