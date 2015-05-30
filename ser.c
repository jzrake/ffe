#include <string.h>
#include "ffe.h"


#if (COW_HDF5)

#include <hdf5.h>

static int read_old_h5type(hid_t src_dset, hid_t dst_type, void *dst_data,
			   const char *type_name);

#endif







int read_write_status(struct ffe_sim *sim, const char *chkpt_name, char mode)
{
  int error = 1;

#if (COW_HDF5)
#define ADD_MEM(nm, tp) H5Tinsert(h5t, #nm, HOFFSET(struct ffe_status, nm), tp)

  if (!H5Fis_hdf5(chkpt_name)) {
    printf("[ser] error: checkpoint '%s' does not exist\n", chkpt_name);
    return 1;
  }

  hid_t h5f = H5Fopen(chkpt_name, H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t h5s = H5Screate(H5S_SCALAR);
  hid_t h5t = H5Tcreate(H5T_COMPOUND, sizeof(struct ffe_status));
  hid_t h5d = -1;

  ADD_MEM(iteration, H5T_NATIVE_INT);
  ADD_MEM(checkpoint_number, H5T_NATIVE_INT);
  ADD_MEM(time_simulation, H5T_NATIVE_DOUBLE);
  ADD_MEM(time_step, H5T_NATIVE_DOUBLE);
  ADD_MEM(time_last_checkpoint, H5T_NATIVE_DOUBLE);
  ADD_MEM(kzps, H5T_NATIVE_DOUBLE);


  if (mode == 'w') {

    if (H5Lexists(h5f, "status", H5P_DEFAULT)) {
      H5Ldelete(h5f, "status", H5P_DEFAULT);
    }

    h5d = H5Dcreate(h5f, "status", h5t, h5s, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(h5d, h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sim->status);

    error = 0;
  }

  else if (mode == 'r') {

    h5d = H5Dopen(h5f, "status", H5P_DEFAULT);

    hid_t src_type = H5Dget_type(h5d);

    if (H5Tequal(src_type, h5t) <= 0) {

      printf("[ser] warning: checkpoint 'status' out of date\n");

      if (read_old_h5type(h5d, h5t, &sim->status, "status")) {
	printf("[ser] error: checkpoint 'status' could not be read\n");
      }
      else {
	error = 0; /* success */
      }

    }
    else {
      H5Dread(h5d, h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sim->status);
      error = 0; /* success */
    }

    H5Tclose(src_type);

  }


  H5Dclose(h5d);
  H5Tclose(h5t);
  H5Sclose(h5s);
  H5Fclose(h5f);

#undef ADD_MEM
#endif
  return error;
}



int read_write_sim(struct ffe_sim *sim, const char *chkpt_name, char mode)
{
  int error = 1;

#if (COW_HDF5)
#define ADD_MEM(nm, tp) H5Tinsert(h5t, #nm, HOFFSET(struct ffe_sim, nm), tp)

  if (!H5Fis_hdf5(chkpt_name)) {
    printf("[ser] error: checkpoint '%s' does not exist\n", chkpt_name);
    return 1;
  }

  hid_t h5t_string_1024 = H5Tcopy(H5T_C_S1); H5Tset_size(h5t_string_1024, 1024);
  hid_t h5f = H5Fopen(chkpt_name, H5F_ACC_RDWR, H5P_DEFAULT);
  hid_t h5s = H5Screate(H5S_SCALAR);
  hid_t h5t = H5Tcreate(H5T_COMPOUND, sizeof(struct ffe_sim));
  hid_t h5d = -1;

  ADD_MEM(Ni, H5T_NATIVE_INT);
  ADD_MEM(Nj, H5T_NATIVE_INT);
  ADD_MEM(Nk, H5T_NATIVE_INT);
  ADD_MEM(cfl_parameter, H5T_NATIVE_DOUBLE);
  ADD_MEM(eps_parameter, H5T_NATIVE_DOUBLE);
  ADD_MEM(kreiss_oliger_mode, H5T_C_S1);
  ADD_MEM(pfeiffer_terms, H5T_C_S1);
  ADD_MEM(time_between_checkpoints, H5T_NATIVE_DOUBLE);
  ADD_MEM(time_final, H5T_NATIVE_DOUBLE);
  ADD_MEM(fractional_helicity, H5T_NATIVE_DOUBLE);
  ADD_MEM(abc_coefficients[0], H5T_NATIVE_DOUBLE);
  ADD_MEM(abc_coefficients[1], H5T_NATIVE_DOUBLE);
  ADD_MEM(abc_coefficients[2], H5T_NATIVE_DOUBLE);
  ADD_MEM(damping_timescale, H5T_NATIVE_DOUBLE);
  ADD_MEM(alpha_squared, H5T_NATIVE_INT);
  ADD_MEM(output_directory, h5t_string_1024);
  ADD_MEM(problem_name, h5t_string_1024);
  ADD_MEM(measure_cadence, H5T_NATIVE_INT);
  ADD_MEM(analyze_cadence, H5T_NATIVE_INT);
  ADD_MEM(num_pspec_bins, H5T_NATIVE_INT);
  ADD_MEM(max_pspec_bin, H5T_NATIVE_INT);

  if (mode == 'w') {

    if (H5Lexists(h5f, "sim", H5P_DEFAULT)) {
      H5Ldelete(h5f, "sim", H5P_DEFAULT);
    }

    h5d = H5Dcreate(h5f, "sim", h5t, h5s, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(h5d, h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, sim);
    error = 0;
  }

  else if (mode == 'r') {

    h5d = H5Dopen(h5f, "sim", H5P_DEFAULT);

    hid_t src_type = H5Dget_type(h5d);

    if (H5Tequal(src_type, h5t) <= 0) {

      printf("[ser] warning: checkpoint 'sim' out of date\n");

      if (read_old_h5type(h5d, h5t, sim, "sim")) {
	printf("[ser] error: checkpoint 'sim' could not be read\n");
      }
      else {
	error = 0; /* success */
      }

    }
    else {
      H5Dread(h5d, h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, sim);
      error = 0; /* success */
    }

    H5Tclose(src_type);

  }

  H5Tclose(h5t_string_1024);
  H5Dclose(h5d);
  H5Tclose(h5t);
  H5Sclose(h5s);
  H5Fclose(h5f);

#undef ADD_MEM
#endif

  return error;
}





#if (COW_HDF5)

static int read_old_h5type(hid_t src_dset, hid_t dst_type, void *dst_data,
			   const char *type_name)
{

  hid_t src_type = H5Dget_type(src_dset);
  int dst_nmembers = H5Tget_nmembers(dst_type);
  int src_nmembers = H5Tget_nmembers(src_type);
  size_t src_size = H5Tget_size(src_type);
  void *src_data = malloc(src_size);

  int error = 0;

  H5Dread(src_dset, src_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, src_data);
      

  for (int n=0; n<dst_nmembers; ++n) {

    char *member_name = H5Tget_member_name(dst_type, n);
    int found = 0;

    //printf("looking for %s in old type\n", member_name);

    for (int m=0; m<src_nmembers; ++m) {

      if (!strcmp(member_name, H5Tget_member_name(src_type, m))) {

	//printf("found %s=%d\n", member_name, m);

	size_t src_member_offset = H5Tget_member_offset(src_type, m);
	size_t dst_member_offset = H5Tget_member_offset(dst_type, n);
	hid_t src_member_type = H5Tget_member_type(src_type, m);
	hid_t dst_member_type = H5Tget_member_type(dst_type, n);

	if (H5Tequal(src_member_type, dst_member_type) <= 0) {

	  printf("[ser] error: incompatible type for member %s:%s\n",
		 type_name, member_name);
	  error += 1;

	}
	else {

	  size_t member_size = H5Tget_size(dst_member_type);
	  memcpy(dst_data + dst_member_offset, src_data + src_member_offset,
		 member_size);

	}

	H5Tclose(src_member_type);
	H5Tclose(dst_member_type);

	found = 1;
	break;
      }
    }

    if (found == 0) {
      printf("[ser] warning: %s:%s not found in source data\n",
	     type_name, member_name);
    }

  }

  free(src_data);
  H5Tclose(src_type);

  return error;
}

#endif
