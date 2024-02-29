&NAM_OASIS NB_TIME_STEPS=<toytimes>,
           DELTA_T=<toydt>,
           GRID_FILENAME='grid_wav.nc' /

&NAM_FCT_SEND CTYPE_FCT='FILES',
              CNAME_FILE='toy_wav.nc',
              VALUE=10 /

&NAM_RECV_FIELDS NB_RECV_FIELDS=2,
                 CRCVFIELDS(1)='TOY_U_01',
                 CRCVFIELDS(2)='TOY_V_01' /

&NAM_SEND_FIELDS NB_SEND_FIELDS=1,
                 CSNDFIELDS(1)='TOY__CHA' /
