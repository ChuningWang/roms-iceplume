      if (typevar == 'uint') then
      i1=2
#  ifdef MPI
      if (WEST_INTER) i1=1
#  endif
      i2=Lmmpi
      j1=1
      j2=Mmmpi
      elseif (typevar == 'u') then
      i1=1
      i2=Lmmpi+1
      j1=0
      j2=Mmmpi+1
      endif
      if (typevar == 'vint') then
      i1=1
      i2=Lmmpi
      j1=2
#  ifdef MPI
      if (SOUTH_INTER) j1=1
#  endif
      j2=Mmmpi
      elseif (typevar == 'v') then
      i1=0
      i2=Lmmpi+1      
      j1=1
      j2=Mmmpi+1
      endif    
      if (typevar == 'rint') then
      i1=1
      i2=Lmmpi
      j1=1
      j2=Mmmpi
      elseif (typevar == 'r') then
      i1=0
      i2=Lmmpi+1
      j1=0
      j2=Mmmpi+1
      endif          