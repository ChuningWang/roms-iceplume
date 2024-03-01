       REAL ::  c1   ,c2    ,c3    ,c4    ,c5    , c6
       REAL :: cb1   ,cb2   ,cb3   ,cb4   ,cb5   ,cbb 
       REAL :: a1    ,a2    ,a3    ,a5    ,nn
       REAL :: ab1   ,ab2   ,ab3   ,ab5   ,nb
       REAL :: sf_d0 ,sf_d1 ,sf_d2 ,sf_d3 ,sf_d4 , sf_d5
       REAL :: sf_n0 ,sf_n1 ,sf_n2
       REAL :: sf_nb0,sf_nb1,sf_nb2                
       REAL :: lim_am0,lim_am1,lim_am2,lim_am3,lim_am4,lim_am5,lim_am6
# ifdef CANUTO_A          /* Canuto & al. A 2001 */
       PARAMETER(c1=5.   ,
     &           c2=0.8  ,
     &           c3=1.968,
     &           c4=1.136,
     &           c5=0.   ,
     &           c6=0.4   )
       PARAMETER(cb1=5.95  ,
     &           cb2=0.6   ,
     &           cb3=1.    ,
     &           cb4=0.    ,
     &           cb5=0.3333,
     &           cbb=0.72   )
# elif defined GibLau_78  /* Gibson Launder 1978 */
       PARAMETER(c1=3.6,
     &           c2=0.8,
     &           c3=1.2,
     &           c4=1.2,
     &           c5=0. ,
     &           c6=0.5 )
       PARAMETER(cb1=3.0,
     &           cb2=0.3333,
     &           cb3=0.333 ,
     &           cb4=0.    ,
     &           cb5=0.3333,
     &           cbb=0.8    )
# elif defined MelYam_82  /* Mellor Yamada  1982 */
       PARAMETER(c1=6.  ,
     &           c2=0.32,  
     &           c3=0.  ,  
     &           c4=0.  ,  
     &           c5=0.  , 
     &           c6=0.   )
       PARAMETER(cb1=3.728,
     &           cb2=0.   , 
     &           cb3=0.   , 
     &           cb4=0.   , 
     &           cb5=0.   , 
     &           cbb=0.6102)
# elif defined KanCla_94  /* Kantha Clayson 1994 */
       PARAMETER(c1=6.  ,
     &           c2=0.32,
     &           c3=0.  ,
     &           c4=0.  ,
     &           c5=0.  ,
     &           c6=0.  )
       PARAMETER(cb1=3.728,
     &           cb2=0.7  ,
     &           cb3=0.7  ,
     &           cb4=0.   ,
     &           cb5=0.2  , 
     &           cbb=0.6102)
# elif defined Luyten_96  /* Luyten & al.   1996 */
       PARAMETER(c1=3. ,
     &           c2=0.8,
     &           c3=2. ,
     &           c4=1.118,
     &           c5=0.   ,
     &           c6=0.5  )
       PARAMETER(cb1=3.,
     &           cb2=0.3333,
     &           cb3=0.3333,
     &           cb4=0.    ,
     &           cb5=0.3333,
     &           cbb=0.8    )
# elif defined CANUTO_B           /* Canuto & al. B 2001 */
       PARAMETER(c1=5.  ,
     &           c2=0.6983,
     &           c3=1.9664,
     &           c4=1.094 ,
     &           c5=0.    ,
     &           c6=0.495  )
       PARAMETER(cb1=5.6,
     &           cb2=0.6,
     &           cb3=1. ,
     &           cb4=0. ,
     &           cb5=0.3333,
     &           cbb=0.477)
# elif defined Cheng_02           /* Cheng 2002 */
       PARAMETER(c1=5.,
     &           c2=0.7983,
     &           c3=1.968 ,
     &           c4=1.136 ,
     &           c5=0.    ,
     &           c6=0.5)
       PARAMETER(cb1=5.52,
     &           cb2=0.2134,
     &           cb3=0.3570,
     &           cb4=0.,
     &           cb5=0.3333,
     &           cbb=0.82)
# else                   /* Canuto & al. A 2001 */
       PARAMETER(c1=5.   ,
     &           c2=0.8  ,
     &           c3=1.968,
     &           c4=1.136,
     &           c5=0.   ,
     &           c6=0.4   )
       PARAMETER(cb1=5.95  ,
     &           cb2=0.6   ,
     &           cb3=1.    ,
     &           cb4=0.    ,
     &           cb5=0.3333,
     &           cbb=0.72   )     
# endif     
       PARAMETER(  a1 = 0.66666666667 - 0.5*c2 )
       PARAMETER(  a2 = 1.            - 0.5*c3 )
       PARAMETER(  a3 = 1.            - 0.5*c4 )            
       PARAMETER(  a5 = 0.5           - 0.5*c6 )
       PARAMETER( ab1 = 1. - cb2               )
       PARAMETER( ab2 = 1. - cb3               )
       PARAMETER( ab3 = 2.*(1.-cb4)            )   
       PARAMETER( ab5 = 2.*cbb*(1.-cb5)        )
       PARAMETER( nn  = 0.5*c1                 )     
       PARAMETER( nb  = cb1                    )       
       PARAMETER( sf_d0 = 36.0*nn*nn*nn*nb*nb                           )
       PARAMETER( sf_d1 = 84.0*a5*ab3*nn*nn*nb+36.0*ab5*nn*nn*nn*nb     )
       PARAMETER( sf_d2 = 9.0*(ab2*ab2-ab1*ab1)*nn*nn*nn 
     &                  - 12.0*(a2*a2-3.*a3*a3)*nn*nb*nb)
       PARAMETER( sf_d3 = 12.0*a5*ab3*(a2*ab1-3.0*a3*ab2)* nn      
     &                    + 12.0*a5*ab3*(    a3*a3-a2*a2)* nb        
     &                    + 12.0*   ab5*(3.0*a3*a3-a2*a2)*nn*nb          )
       PARAMETER( sf_d4 = 48.0*a5*a5*ab3*ab3*nn + 36.0*a5*ab3*ab5*nn*nn )
       PARAMETER( sf_d5 = 3.0*(a2*a2-3.0*a3*a3)
     &                       *(ab1*ab1-ab2*ab2)*nn    )
       PARAMETER( sf_n0  = 36.0*a1*nn*nn*nb*nb )
       PARAMETER( sf_n1  = - 12.0*a5*ab3*(ab1+ab2)*nn*nn           
     &                    + 8.0*a5*ab3*(6.0*a1-a2-3.0*a3)*nn*nb     
     &                    + 36.0*a1*ab5*nn*nn*Nb )
       PARAMETER( sf_n2  = 9.0*a1*(ab2*ab2-ab1*ab1)*nn*nn ) 
       PARAMETER( sf_nb0 = 12.0*ab3*nn*nn*nn*nb  )
       PARAMETER( sf_nb1 = 12.0*a5*ab3*ab3*nn*nn )
       PARAMETER( sf_nb2 = 9.0*a1*ab3*(ab1-ab2)*nn*nn + ( 6.0*a1*(a2-3.0*a3) 
     &                               - 4.0*(a2*a2-3.0*a3*a3) )*ab3 * nn * nb)
       PARAMETER( lim_am0 = sf_d0*sf_n0               )
       PARAMETER( lim_am1 = sf_d0*sf_n1 + sf_d1*sf_n0 )
       PARAMETER( lim_am2 = sf_d1*sf_n1 + sf_d4*sf_n0 )
       PARAMETER( lim_am3 = sf_d4*sf_n1               )
       PARAMETER( lim_am4 = sf_d2*sf_n0               )
       PARAMETER( lim_am5 = sf_d2*sf_n1+sf_d3*sf_n0   )
       PARAMETER( lim_am6 = sf_d3*sf_n1               )      
     
     
