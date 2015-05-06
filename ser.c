#include "ffe.h"

#if (COW_HDF5)
#include <hdf5.h>
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

    hid_t type = H5Dget_type(h5d);

    if (H5Tequal(type, h5t) <= 0) {
      printf("[ser] error: checkpoint format 'status' is not up-to-date\n");
    }
    else {
      H5Dread(h5d, h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sim->status);
      error = 0;
    }

    H5Tclose(type);
 
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
  ADD_MEM(alpha_squared, H5T_NATIVE_INT);
  ADD_MEM(output_directory, h5t_string_1024);
  ADD_MEM(problem_name, h5t_string_1024);

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

    hid_t type = H5Dget_type(h5d);

    if (H5Tequal(type, h5t) <= 0) {
      printf("[ser] error: checkpoint format 'sim' is not up-to-date\n");
    }
    else {
      H5Dread(h5d, h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, sim);
      error = 0;
    }

    H5Tclose(type);

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
