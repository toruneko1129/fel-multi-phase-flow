NAME	=	a.out

FFLAGS=-Kfast,parallel,openmp
FC = mpifrtpx

SRCDIR	=	./srcs
FCSRCS =            \
bnd_avzero.f			 \
bnd_contact_angle.f \
bnd_dirichlet.f                  \
bnd_comm.f			 \
bnd_neumann.f			 \
bnd_periodic.f			 \
bndu.f				 \
cal_advu.f			 \
cal_advu_weno5.f		 \
cal_arith_coef_vis.f		 \
cal_arith_mean.f		 \
cal_arith_tau.f			 \
cal_av.f			 \
cal_center.f                     \
cal_coef_prs.f			 \
cal_coef_solp_fft.f		 \
cal_div_tensor.f		 \
cal_grad_p2a.f			 \
cal_grad_p2uvw.f		 \
cal_grad_prs.f			 \
cal_harm_coef_vis.f		 \
cal_harm_mean.f			 \
cal_harm_tau.f			 \
cal_srcu.f			 \
cal_vel.f                        \
caldiv.f                         \
caldt.f				 \
calfst.f			 \
calphat.f			 \
calrkap.f			 \
calq.f                           \
calrn.f				 \
calsij.f			 \
caltau.f			 \
calvorx.f                        \
corunp_explicit.f		 \
cpy.f				 \
datain.f			 \
dataou.f			 \
fourn.f                          \
init.f				 \
init_q.f                         \
mk_all.f			 \
mk_center.f                      \
mk_slice_k.f			 \
mk_slice_j.f                     \
mk_vel.f                         \
mkvtk_phi.f			 \
mkvtk_phil.f                      \
mkvtk_q.f			 \
mkvtk_p.f                        \
set_ransu192.f                   \
solp_fft_tdma1.f		 \
solp_fft_tdma2.f		 \
solp_fft_tdma3.f		 \
solp_fft_tdma4.f		 \
solphi_mthinc1.f		 \
solphi_mthinc2.f		 \
solphi_mthinc3.f		 \
solu_sor4.f			 \
summation.f                      \
sum_fst.f                        \
trans_r2w.f			 \
trans_w2r.f                      \

F90SRCS =           \
init_phi_legendre.f90 \
main.f90            \

SRCS	=	$(addprefix $(SRCDIR)/, $(F90SRCS)) \
			$(addprefix $(SRCDIR)/, $(FCSRCS))

OBJDIR	=	./objs
OBJS	=	$(F90SRCS:%.f90=$(OBJDIR)/%.o) \
			$(FCSRCS:%.f=$(OBJDIR)/%.o)

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
		mkdir -p $(OBJDIR)
		$(FC) -c $(FFLAGS) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.f
		mkdir -p $(OBJDIR)
		$(FC) -c $(FFLAGS) $< -o $@

$(NAME): $(OBJS)
	$(FC) $(OBJS) $(FFLAGS) -o $(NAME)

all: $(NAME)

clean:
	rm -rf $(OBJDIR)

fclean: clean
	rm -rf $(NAME)

re: fclean all

.SUFFIXES: .o .f .f90