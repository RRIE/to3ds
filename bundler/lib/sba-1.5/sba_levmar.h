typedef struct wrap_motstr_data_ {
  void   (*proj)(int j, int i, double *aj, double *bi, double *xij, void *adata); // Q
  void (*projac)(int j, int i, double *aj, double *bi, double *Aij, double *Bij, void *adata); // dQ/da, dQ/db
  int cnp, pnp, mnp; /* parameter numbers */
  void *adata;
};
