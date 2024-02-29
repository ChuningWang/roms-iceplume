&NAM_OASIS NB_TIME_STEPS=<toytimes>,
           DELTA_T=<toydt>,
           GRID_FILENAME='grid_atm.nc' /

&NAM_FCT_SEND CTYPE_FCT='FILES',
              CNAME_FILE='toy_atm.nc',
              VALUE=10 /

&NAM_RECV_FIELDS NB_RECV_FIELDS=1,
                 CRCVFIELDS(1)='TOY__CHA' /

&NAM_SEND_FIELDS NB_SEND_FIELDS=2,
                 CSNDFIELDS(1)='TOY_U_01',
                 CSNDFIELDS(2)='TOY_V_01' /
