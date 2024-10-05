cc
      dfm_1= dfm5*2.0d0 - dfm3*7.0d0 + dfm1*1.1d1
      dfm_2=-dfm3       + dfm1*5.0d0 + dfp1*2.0d0
      dfm_3= dfm1*2.0d0 + dfp1*5.0d0 - dfp3

      dfp_1= dfp5*2.0d0 - dfp3*7.0d0 + dfp1*1.1d1
      dfp_2=-dfp3       + dfp1*5.0d0 + dfm1*2.0d0
      dfp_3= dfp1*2.0d0 + dfm1*5.0d0 - dfm3

      sm_1=(dfm5       - dfm3*2.0d0 + dfm1      )**2*1.3d1
     1    +(dfm5       - dfm3*4.0d0 + dfm1*3.0d0)**2*3.0d0
      sm_2=(dfm3       - dfm1*2.0d0 + dfp1      )**2*1.3d1
     1    +(dfm3                    - dfp1      )**2*3.0d0
      sm_3=(dfm1       - dfp1*2.0d0 + dfp3      )**2*1.3d1
     1    +(dfm1*3.0d0 - dfp1*4.0d0 + dfp3      )**2*3.0d0

      sp_1=(dfp5       - dfp3*2.0d0 + dfp1      )**2*1.3d1
     1    +(dfp5       - dfp3*4.0d0 + dfp1*3.0d0)**2*3.0d0
      sp_2=(dfp3       - dfp1*2.0d0 + dfm1      )**2*1.3d1
     1    +(dfp3                    - dfm1      )**2*3.0d0
      sp_3=(dfp1       - dfm1*2.0d0 + dfm3      )**2*1.3d1
     1    +(dfp1*3.0d0 - dfm1*4.0d0 + dfm3      )**2*3.0d0

      am_1=sm_2*sm_3
      am_2=sm_1*sm_3
      am_3=sm_1*sm_2

      ap_1=sp_2*sp_3
      ap_2=sp_1*sp_3
      ap_3=sp_1*sp_2

ccc   am_1:(sm_2**2)*(sm_3**2)
ccc   am_2:6*(sm_1**2)*(sm_3**2)
ccc   am_3:3*(sm_1**2)*(sm_2**2)

      am_1=am_1**2		
      am_2=am_2**2
      am_3=am_3**2

      ap_1=ap_1**2
      ap_2=ap_2**2
      ap_3=ap_3**2
      
      am_2=am_2*6.0d0
      am_3=am_3*3.0d0

      ap_2=ap_2*6.0d0
      ap_3=ap_3*3.0d0

      deninv_am=1.0d0/(am_1+am_2+am_3+verysmall)
      deninv_ap=1.0d0/(ap_1+ap_2+ap_3+verysmall)

      wgm_1=am_1*deninv_am
      wgm_3=am_3*deninv_am

      wgp_1=ap_1*deninv_ap
      wgp_3=ap_3*deninv_ap

      wgm_2=1.0d0-wgm_1-wgm_3
      wgp_2=1.0d0-wgp_1-wgp_3

      velab=abs(vel00)

      velm1=vel00+velab
      velp1=vel00-velab

      adv_weno5=(
     1    velm1*(wgm_1*dfm_1 + wgm_2*dfm_2 + wgm_3*dfm_3)
     2   +velp1*(wgp_1*dfp_1 + wgp_2*dfp_2 + wgp_3*dfp_3)
     & )*d1_12
cc
